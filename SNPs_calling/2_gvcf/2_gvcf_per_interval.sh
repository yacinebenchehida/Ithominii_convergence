#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=SNP_CAL
#SBATCH --array=0-59

################################
# 0 - DEFINE WORKING DIRECTORY #
################################
path_bam="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Results/1_sorted_dedup_bam/$1"
intervals="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$2/*0${SLURM_ARRAY_TASK_ID}-scattered.interval_list"
path_result="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/2_gvcf/Results/$1"
ref_genome="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$2"
fasta_sequence=$(ls $ref_genome|grep -E "*.fa$|*.fasta$")

##################################
# 0bis - Load necessar softwares #
##################################
module load  GATK/4.3.0.0-GCCcore-11.3.0-Java-11

mkdir -p  $path_result

##########################################################
# 1 Run haplotype called and get gvcfs for each interval #
##########################################################
gatk --java-options "-Xmx4g" HaplotypeCaller \
-R $ref_genome/$fasta_sequence \
-I $path_bam/*bam \
-O $path_result/"$1"_${SLURM_ARRAY_TASK_ID}.g.vcf.gz \
--intervals $intervals \
--output-mode EMIT_ALL_CONFIDENT_SITES \
-ERC GVCF \
--dont-use-soft-clipped-bases \
--native-pair-hmm-threads 4
