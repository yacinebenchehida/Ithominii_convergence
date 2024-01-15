#! /bin/bash

#################
# Set variables #
#################
INPUT="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/Results/$1/braker_combined_tsebra.gtf"
START=$2
END=$3
SCAFFOLD=$4

##################################
# Extract annotation around peak #
##################################
cat $INPUT |grep -E $SCAFFOLD |awk -v start="$START" -v end="$END" '($4 >= start) && ($5 < end)'|grep -P "\tgene\t|exon|intron|codon"|perl -pe 's/(.+)transcript_id(.+)(g_\d+).+/$1$3/g'|awk '{print $9"\t"$3"\t"$4"\t"$5"\t"$7}' > Gene_position.txt

################################################################################
# Clean annotation to plot by making sure that the genes selected are complete #
################################################################################
module load R/4.2.1-foss-2022a
Rscript Complete_genes_selection.R $1
