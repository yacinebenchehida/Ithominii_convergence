#!/bin/bash
#Author: Yacine Ben Chehida

#SBATCH --time=1-06:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --account=BIOL-SPECGEN-2018

#############
# Set paths #
#############
RESULTS="/users/ybc502/scratch/Conv_Evol/Analyses_bio/Ref_genome_Cleaning/Blast_per_window/Results/$2"
INPUT="/users/ybc502/scratch/Conv_Evol/Analyses_bio/Ref_genome_Cleaning/Blast_per_window/Inputs"
SCRIPT="/users/ybc502/scratch/Conv_Evol/Analyses_bio/Ref_genome_Cleaning/Blast_per_window/Script"
USER_NAME="ybc502"
NTDB="/mnt/lustre/groups/biol-specgen-2018/yacine/data/genbank/ntdb"

##################
# Load biopython #
##################
module load bio/Biopython/1.79-foss-2021a

#######################################
# Calculate the size of each scaffold #
#######################################
mkdir -p $INPUT/$2
python -W ignore scaffold_size.py $1 > $INPUT/$2/Scaffolds_size.txt
echo Scaffold sizes calculated

#########################################################
# Cut the scaffolds into windows of user specified size #
#########################################################
## 1) Create a separate folder for each scaffold
cat $INPUT/$2/Scaffolds_size.txt|awk '{print $1}' | while read line; do mkdir -p $INPUT/$2/$line ; done
echo A folder for each scaffold created

## 2) Create a fasta file for each scaffold
python -W ignore Fasta_per_scaffold.py $1 $INPUT/$2
echo A fasta file for each scaffold created

## 3) Chop each scaffold into a constant size windows
cat $INPUT/$2/Scaffolds_size.txt|while read line; do scaffold=$(echo $line|awk '{print $1}'); size=$(echo $line|awk '{print $2}'); sbatch --job-name=win_fas ./Fasta_windows_chopper.sh $scaffold $3 $size $2; done

## 4) Extract the running job numbers
running_jobs1=$(squeue|grep $USER_NAME| grep win_fas| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')

## 5) Create a window interval file. Necessary for a later step.
cat $INPUT/$2/Scaffolds_size.txt|while read line; do scaffold=$(echo $line|awk '{print $1}'); size=$(echo $line|awk '{print $2}'); ./windows_file_creator.sh $scaffold $3 $size $2 $INPUT; done
echo faster interval file created for all scaffolds

##########################################
# Blasting each window for each scaffold #
##########################################
cd $INPUT/$2

for i in $(cat Scaffolds_size.txt|awk '{print $1}') 
	do cd $i
	echo $i
	sleep 2m
	mkdir -p $RESULTS/$i
	echo -e qseqid"\t"staxids"\t"sseqid"\t"pident"\t"length"\t"mismatch"\t"gapopen"\t"qstart"\t"qend"\t"sstart"\t"send"\t"evalue"\t"bitscore"\t"TAXA"\t"window_start"\t"window_end > $RESULTS/$i/Blast_hits.txt
	NB_ARRAY=$(cat *_windows_position.txt | wc -l)
	sbatch --array=1-$NB_ARRAY --dependency=aftercorr:$running_jobs1 $SCRIPT/blast_script.sh $i $2 $INPUT $RESULTS $NTDB
	cd ..
done

# --dependency=aftercorr:$running_jobs1
###################
# Combine results #
###################
running_jobs2=$(squeue|grep $USER_NAME| grep blast_sc| awk '{print $1}'|perl -pe 's/(\d+)_(.+)/$1/g'|sort|uniq|perl -pe 's/\n/,/g'|sed 's/,$//g')
echo $running_jobs2
sbatch --dependency=aftercorr:$running_jobs2 $SCRIPT/Combined_results.sh $RESULTS
