#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=Filt_vcfs
#SBATCH --array=0-59

###########################
# 0 - DEFINE USEFUL PATHS #
###########################
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/4_Genotype_gvcfs/Results/$1"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/$1"


################################
# 1 - Load necessary softwares #
################################
module load  BCFtools/1.15-GCC-11.2.0

#############################  
# 2 - SET WORKING DIRECTORY #
#############################
cd  /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/4_Genotype_gvcfs/

#############################   
# 3 - CREATE RESULTS FOLDER #
#############################   
mkdir -p $RESULTS

####################################################################  
# 4 - FILTER VCF USING BCFTOOLS TO KEEP ONLY GOOD QUALITY VARIANTS #
####################################################################  
bcftools filter -e 'FORMAT/DP < 1 |FORMAT/GQ < 5 |QUAL <= 5' --set-GTs . $VCF/"$1"_genotypeGVCF.intervals_${SLURM_ARRAY_TASK_ID}.vcf.gz -O u | bcftools view -U -i 'TYPE=="snp"' -m2 -M2 -v snps -O v| bcftools view -i 'F_MISSING < 0.2'> $RESULTS/"$1"_genotypeGVCF.intervals_${SLURM_ARRAY_TASK_ID}.filters.DP4_GQ19_QUAL10_F_M25.snps.vcf
