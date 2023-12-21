#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=35G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=map_sort_RG_dup_merge

################################
# 0 - DEFINE WORKING DIRECTORY #
################################
path_fastq_files="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/Fastq_files/$3/$1"
path_results="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Results"
ref_genome="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$2"
fasta_sequence=$(ls $ref_genome|grep -E "*.fa$|*.fasta$")

##################################
# 0bis - Load necessar softwares #
##################################
module load BWA/0.7.17-GCCcore-11.2.0 # Load bwa
module load picard/2.25.5-Java-13 # Load picard
module load SAMtools/1.16.1-GCC-11.3.0 # Load Samtools
module load Qualimap/2.2.1-foss-2019b-R-3.6.2

###########################################################
# 1 - Align reads + transform sam to bam + sort bam files #
###########################################################
mkdir -p $path_results/0_sorted_bam/$1 # Create working directory for sorted bam files
bwa mem -t 12 $ref_genome/$fasta_sequence $path_fastq_files/*_R1*.fastq.gz $path_fastq_files/*R2*.fastq.gz | samtools view -bSh |samtools sort -@ 12 -o $path_results/0_sorted_bam/$1/"$1"_sorted.bam # Use bwa, samtools view and samtools sort

#############################
# 2 - Picard ADD read group #
#############################
mkdir -p $path_results/1_sorted_dedup_bam/$1  # Create working directory for merged bam file without duplicates
mkdir -p /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Scripts/TMP/$1

RGID=$(zcat $path_fastq_files/*R1*.fastq.gz | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4) # Extract RGID from read header
RGLB="Solexa"
RGSM="$1"
RGPU=$RGID.$LB
RGPL="Illumina" 

java -jar -Xmx50g $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
 I=$path_results/0_sorted_bam/$1/"$1"_sorted.bam \
 O=$path_results/0_sorted_bam/$1/"$1"_sorted_RG.bam \
 RGID=$RGID \
 RGLB=$RGLB \
 RGPL=$RGPL \
 RGPU=$RGPU \
 RGSM=$RGSM \
 SORT_ORDER=coordinate \
 CREATE_INDEX=TRUE

#############################
# 3 - Picard MarkDuplicates #
#############################
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=/users/ybc502/scratch/Conv_Evol/TMP/$1 \
I=$path_results/0_sorted_bam/$1/"$1"_sorted_RG.bam \
O=$path_results/1_sorted_dedup_bam/$1/"$1"_sorted_dedup.bam \
M=$path_results/1_sorted_dedup_bam/$1/"$1"_marked_dedup_metrics.txt

rm $path_results/0_sorted_bam/$1/*.bam # Remove heavy sorted bam with read groups

#######################################
# 4 - generating  flagstat then index #
#######################################
samtools flagstat -@12 $path_results/1_sorted_dedup_bam/$1/"$1"_sorted_dedup.bam -O tsv > $path_results/1_sorted_dedup_bam/$1/"$1".flagstat
samtools index -@12 $path_results/1_sorted_dedup_bam/$1/"$1"_sorted_dedup.bam

###########################################                              
# 5 - check mapping quality with qualimap #
###########################################
mkdir -p $path_results/2_Qualimap/"$1"

qualimap bamqc --java-mem-size=20G \
-bam $path_results/1_sorted_dedup_bam/$1/"$1"_sorted_dedup.bam \
-c \
-gd HUMAN \
-nt 12 \
-outdir $path_results/2_Qualimap/"$1" \
-outformat HTML 
