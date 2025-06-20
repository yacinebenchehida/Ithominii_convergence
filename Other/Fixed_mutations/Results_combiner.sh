#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-00:25:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G 
#SBATCH --account=BIOL-SPECGEN-2018 

################
# Define paths #
################
PATH_RESULTS="/mnt/lustre/groups/biol-specgen-2018/yacine/Fixed_mutations_sliding_windows/Results/$1"
PATH_SCRIPT=$(pwd)

###########################
# Go to working directory #
###########################
cd $PATH_RESULTS

###################################     
# Combine fixed mutations results #
###################################   
echo -e Contig"\t"Start"\t"Stop"\t"Fixed"\t"Het_Hom > "$1"_diff_het_hom.txt
ls *_fixed_mutations.txt|sort -V| xargs cat |grep -v CHROM >> "$1"_diff_het_hom.txt

############################
# Remove all the tmp files #
############################
rm *fixed_mutations.txt

#####################################################
# Plot fixed and het/hom mutations along the genome #
#####################################################
module load lang/R/4.2.1-foss-2022a
Rscript --vanilla $PATH_SCRIPT/Plot_fixed_mutation_along_genome.R $1_diff_het_hom.txt $1 $2 $PATH_RESULTS
