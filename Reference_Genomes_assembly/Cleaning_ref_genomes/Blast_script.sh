#!/bin/bash
#Author: Yacine Ben Chehida

#SBATCH --time=0-03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --account=BIOL-SPECGEN-2018


#############
# Set paths #
#############
NTDB="$5"
SCAFF=$1
RESULTS="$4/$SCAFF"
INPUT="$3/$2/$SCAFF"

###########################
# Load the blast+ library #
###########################
module load bio/Biopython/1.79-foss-2021a

################# 
# Select window #
#################
WINDOW_START=$(cat $INPUT/*windows_position.txt | sed -n "${SLURM_ARRAY_TASK_ID}"p | awk '{print $1}')
WINDOW_END=$(cat $INPUT/*windows_position.txt | sed -n "${SLURM_ARRAY_TASK_ID}"p | awk '{print $2}')

################################
# Blastn the selected scaffold #
################################
module load bio/BLAST+/2.13.0-gompi-2022a

cd $NTDB
blastn -query $INPUT/"$SCAFF"_"$WINDOW_START"_"$WINDOW_END"/"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".fas -db ntdb -outfmt '6 qseqid staxids sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1 -max_hsps 1  > $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".txt 

#################################################################
# Check taxonomy and append the information to the blast result #
#################################################################
# extract taxid
TAXID=$(cat $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".txt |awk '{print $2}')

# get taxonomy with edirect
/users/ybc502/edirect/efetch -db taxonomy -id $TAXID -format xml | xtract -pattern Taxon -first TaxId -element Taxon -block "*/Taxon" -unless Rank -equals "no rank" -tab " " -sep ":" -element Rank,ScientificName > $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".taxo

# Check if it's a lepidoptera, a bacteria or something else 
TAXON=$(cat $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".taxo |perl -pe 's/(.+)(\s)order:(\w+) (.+)/$3/g')

if [[ $TAXON == "Lepidoptera" ]]; then
        paste <(cat $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".txt) <(echo -e $TAXON"\t"$WINDOW_START"\t"$WINDOW_END)  >>  $RESULTS/Blast_hits.txt
else 
        TAXON=$(cat $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".taxo |perl -pe 's/(.+)(\s)superkingdom:(\w+) (.+)/$3/g')
                if [[ $TAXON == "Bacteria" ]]; then
                        echo $TAXON
                        paste <(cat $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".txt) <(echo -e $TAXON"\t"$WINDOW_START"\t"$WINDOW_END) >> $RESULTS/Blast_hits.txt
                else
                        TAXON=$(cat $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".taxo |perl -pe 's/(.+)(\s)class:(\w+) (.+)/$3/g')
                        echo $TAXON
                        paste <(cat $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".txt) <(echo -e $TAXON"\t"$WINDOW_START"\t"$WINDOW_END) >> $RESULTS/Blast_hits.txt 
                fi
fi

###################################
# Clean (erase) useless tmp files #
###################################
rm $RESULTS/Best_blast_hit_"$SCAFF"_"$WINDOW_START"_"$WINDOW_END".*
#rm -r $INPUT/"$SCAFF"_"$WINDOW_START"_"$WINDOW_END"
#test -f $INPUT/"$SCAFF".fasta && rm $INPUT/"$SCAFF".fasta
