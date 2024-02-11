#!/bin/bash

#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=bam2fas
#SBATCH --time=0-01:00:00
#SBATCH --array=1-40

##########################
# load necessary modules #
##########################
module load R/4.2.1-foss-2022a
module load Biopython/1.81-foss-2022b
module load HTSlib/1.15.1-GCC-11.3.0
module load Java/13.0.2

##########################
# Set useful directories #
##########################
BAM="/mnt/scratch/projects/biol-specgen-2018/edd/phylogeny/aligned_methods/1_mapping/Results/1_sorted_dedup_bam/"
SAMPLES=$(ls $BAM|grep Sample|grep -v -E "Sample_15-2002-658|Sample_11-2002-3511|Sample_2-2002-1553|Sample_2-2002_1942_M_ma_phasiana|Sample_39-2022-0522|Sample_61-2002-3158")
SAMPLE=$(echo "$SAMPLES" | awk -v idx=${SLURM_ARRAY_TASK_ID} 'NR == idx')
echo $SAMPLE
INPUT="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/ASTRAL/Inputs"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/ASTRAL/Results"
TRIMAL="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/trimal/source/trimal"
ANGSD="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/angsd/angsd"
ASTRAL="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/ASTRAL/Astral/astral.5.7.3.jar"

###############################################
# For each scaffold convert fasta to bam file #
###############################################
# list of bam path (with the bam files to analyse)
BAM_FILE=$(readlink -f  $BAM/$SAMPLE/*bam )
echo bam file is $BAM_FILE

# Transform bam to fasta
SCAFFOLD="scaffold_{1..25}"

for scaffold in $(eval echo "$SCAFFOLD"); do
	echo $scaffold
	mkdir -p $RESULTS/$scaffold
	$ANGSD -i  $BAM_FILE -doFasta 3 -out $RESULTS/$scaffold/"$SAMPLE"_"$scaffold" -doCounts 1 -basesPerLine 1000000000 -r $scaffold
	rm $RESULTS/$scaffold/$SAMPLE*arg
done
