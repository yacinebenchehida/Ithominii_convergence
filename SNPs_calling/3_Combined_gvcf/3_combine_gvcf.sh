#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=Comb_gvcf
#SBATCH --array=0-59
#SBATCH --mail-type=END,FAIL

###########################
# 0 - DEFINE USEFUL PATHS #
###########################
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/3_Combined_gvcf/Results/$1"
METADATA="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/3_Combined_gvcf/Scripts"
REF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$1"
FASTA_REF=$(ls $REF|grep -E "*.fa$|*.fasta$")

################################
# 1 - Load necessary softwares #
################################
module load  GATK/4.3.0.0-GCCcore-11.3.0-Java-11

#############################  
# 2 - SET WORKING DIRECTORY #
#############################
cd /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/2_*/Results/

#############################   
# 3 - CREATE RESULTS FOLDER #
#############################   
mkdir -p $RESULTS

################################################   
# 4 - CREATE AUTOMATICALLY THE VARIANT COMMAND #
################################################
variants=$(for i in $(cat $METADATA/ref_genome.txt|grep -v -E  "H_s|H_a"|grep $1|awk '{print $1}'); do cd $i; echo "--variant $i/"$(ls "$i"_${SLURM_ARRAY_TASK_ID}.g.vcf.gz); cd ..; done | perl -pe 's/\n/ /g')

###############################################################  
# 5 - COMBINE GVCFS SEPARATELY PER INTERVALS FOR EACH SPECIES #
###############################################################  
gatk --java-options -Xmx15g CombineGVCFs -R $REF/$FASTA_REF $variants -O $RESULTS/$1_intervals_merged_${SLURM_ARRAY_TASK_ID} -intervals $REF/*0${SLURM_ARRAY_TASK_ID}-scattered.interval_list 
