#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=PCA

###########################
# 1 - DEFINE USEFUL PATHS #
###########################
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/9_Phasing/Results/$1/*vcf.gz"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/PCA/Results"
PHENO="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/PCA/Inputs/$1/grouping_subspecies.txt"
PLINK="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/plink_linux_x86_64_20231018/plink"

########################   
# 2 - LOAD PLINK and R #
######################## 
module load R/4.2.1-foss-2022a
module load BCFtools/1.15.1-GCC-11.3.0
mkdir -p $RESULTS/$1

######################################
# 3 - CREATE A LINKAGE FREE DATA SET #
######################################
$PLINK --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 10 0.1 --out $RESULTS/$1/$1 --geno 0 # set a window of 100 Kb. The second argument, 10 is our window step size - meaning we move 10 bp each time we calculate linkage. Finally, we set an r2 threshold - i.e. the threshold of linkage we are willing to tolerate. Here we prune any variables that show an r2 of greater than 0.05.

#######################
# 4 - PERFORM THE PCA #
#######################
# LD prunning
$PLINK --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --extract $RESULTS/$1/*prune.in --out $RESULTS/$1/$1 --geno 0
# No LD pruning
$PLINK --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --out $RESULTS/$1/$1 -geno 0

####################################
# 5 - ERASE BIG INTERMEDIATE FILES #
####################################
rm  $RESULTS/$1/*log $RESULTS/$1/*sex

########################################################
# 6 - Plot PCA scatter plot and pdf and an html widget #
########################################################
Rscript --vanilla ./PCA_plot.R $RESULTS/$1/*vec $RESULTS/$1/*val $PHENO $RESULTS/$1 $1

# Usage: sbatch ./PCA.sh Melinaea_menophilus
