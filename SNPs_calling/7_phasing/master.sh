#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=master

#############
# SET PATHS #
#############
VCF="/mnt/scratch/projects/biol-specgen-2018/edd/H_anastasia/6_combine_intervals/Results/"$1".basic.filters.DP5_GQ10_Q10.snps.vcf.gz"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/9_Phasing/Results/$1"
REF_PATH="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$1"
REF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$1"
FASTA_REF=$(ls $REF|grep -E "*.fa$|*.fasta$|*.fna$")

################
# Load modules #
################
module load HTSlib/1.15.1-GCC-11.3.0
module load BCFtools/1.15.1-GCC-11.3.0
module load Biopython/1.81-foss-2022b

############################
# Create working directory #
############################ 
mkdir -p $RESULTS

#####################
# GET SCAFFOLD NAMES #
#####################
python scaffold_name.py $REF/$FASTA_REF > order_scaff.txt


#################################
# Get total number of scaffolds #
#################################
array_number=$(cat order_scaff.txt|wc -l)

##################################
# Run shapeit for each scaffolds #
##################################
sbatch --array=1-$array_number ./shapeit.sh $1 $VCF

########################################### 
# Merge all phased vcf into a single file #
###########################################
running_jobs1=$(squeue|grep ybc502| grep phase| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
sbatch --dependency=aftercorr:$running_jobs1 ./combine.sh $array_number
