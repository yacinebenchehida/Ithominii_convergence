#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=Gen_gvcf
#SBATCH --array=1-59

###########################
# 0 - DEFINE USEFUL PATHS #
###########################
GVCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/3_Combined_gvcf/Results/$1"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/4_Genotype_gvcfs/Results/$1"
REF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$1"
FASTA_REF=$(ls $REF|grep -E "*.fa$|*.fasta$")

################################
# 1 - Load necessary softwares #
################################
module load  GATK/4.3.0.0-GCCcore-11.3.0-Java-11

#############################  
# 2 - SET WORKING DIRECTORY #
#############################
cd /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/3_Combined_gvcf/Results/

#############################   
# 3 - CREATE RESULTS FOLDER #
#############################   
mkdir -p $RESULTS

###############################################################  
# 4 - GENOTYPE GVCFS SEPARATELY PER INTERVALS FOR EACH SPECIES #
###############################################################  
gatk --java-options -Xmx15g GenotypeGVCFs -R $REF/$FASTA_REF -V $GVCF/"$1"_intervals_${SLURM_ARRAY_TASK_ID} -O $RESULTS/"$1"_genotypeGVCF.intervals_${SLURM_ARRAY_TASK_ID}.vcf.gz -intervals $REF/*0${SLURM_ARRAY_TASK_ID}-scattered.interval_list 
