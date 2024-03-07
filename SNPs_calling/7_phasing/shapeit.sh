#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=phase

###########################
# 1 - DEFINE USEFUL PATHS #
###########################
VCF="/mnt/scratch/projects/biol-specgen-2018/edd/H_anastasia/5_filtering/Results/"$1"
#VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/$1/"$1".DP4_GQ5_QUAL5_F20.snps.vcf.gz"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/9_Phasing/Results/$1"
SHAPEIT="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/shapeit4/bin/shapeit4.2"
REF_PATH="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$1"
module load HTSlib/1.15.1-GCC-11.3.0
module load module load BCFtools/1.15.1-GCC-11.3.0
module load Biopython/1.81-foss-2022b

#######################################
# 2 - EXTRACT RIGHT INTERVAL POSITION #
#######################################
SCAFFOLD=$(sed -n "${SLURM_ARRAY_TASK_ID}"p order_scaff.txt)

####################################  
# 3 - PHASE THE DATA USING SHAPE IT #
#####################################
$SHAPEIT --input $VCF --region $SCAFFOLD --output $RESULTS/"$1"_genotypeGVCF.intervals_${SLURM_ARRAY_TASK_ID}.filters.DP4_GQ5_QUAL5_FM20.snps.phased.vcf --thread 8

#################################  
# 4 - bgzip and tabix vcf files #
#################################
bgzip $RESULTS/"$1"_genotypeGVCF.intervals_${SLURM_ARRAY_TASK_ID}.filters.DP4_GQ5_QUAL5_FM20.snps.phased.vcf
tabix $RESULTS/"$1"_genotypeGVCF.intervals_${SLURM_ARRAY_TASK_ID}.filters.DP4_GQ5_QUAL5_FM20.snps.phased.vcf.gz
