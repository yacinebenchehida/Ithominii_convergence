#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=5-10:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=ST_BR
#SBATCH --partition=week

###########################
# 0 - Load useful modules #
###########################
module load Biopython/1.79-foss-2022a
module load STAR/2.7.10b-GCC-11.3.0
module load SAMtools/1.17-GCC-12.2.0
module load BRAKER/2.1.6-foss-2022a
module load RepeatModeler/2.0.4-foss-2022a

#############################
# 1 - Set reference genomes #
#############################
DATA="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/Data/$1"
FASTA_REF=$(ls $DATA|grep -E "*.fa$|*.fasta$")
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/Results/$1"

echo $DATA/$FASTA_REF

#####################
# 2 - Repeat masker #
#####################
cd $DATA
BuildDatabase -name $DATA/lepidoptera $DATA/$FASTA_REF 
RepeatModeler -database lepidoptera -threads 8 -LTRStruct 
echo repeatmodeler ran
RepeatMasker -s -xsmall -pa 8 -gff -lib lepidoptera-families.fa $DATA/$FASTA_REF
REF_GENOME="$DATA/*.masked"
echo $REF_GENOME
echo RepeatMasker ran

####################################
# 3 - Build genome index with STAR #
####################################
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $RESULTS/index --genomeFastaFiles $REF_GENOME  --genomeSAindexNbases 13
echo STAR INDEX BUILT

###########################################
# 4 - Align reads to reference  with STAR #
########################################### 
STAR --genomeDir $RESULTS/index --runThreadN 12 --readFilesIn $DATA/RNA_*_R1.fastq $DATA/RNA_*_R2.fastq --alignIntronMax 500000 --outFileNamePrefix $RESULTS/$1 --twopassMode Basic
echo STAR ALIGNMENT PERFORMED

##########################
# 5 - Convert SAM to BAM #
##########################
samtools view -bS "$RESULTS"/*.out.sam|samtools sort -@ 12 -o "$RESULTS"/"$1".bam
echo SAM CONVERT TO BAM FILE

##################
# 6 - Run BRAKER #
##################
REF_GENOME=$(ls $DATA/*.fasta.masked)
echo $REF_GENOME
cd $RESULTS
# First braker run (based on RNA-seq data)
braker.pl --genome=$REF_GENOME --bam="$RESULTS"/"$1".bam --softmasking --cores=12  --AUGUSTUS_CONFIG_PATH=/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/Augustus/config --round 50
# Second braker run (based on protein data and using hintsfile from previous run; i.e. based on step 1 training)
braker.pl --genome=$REF_GENOME  --species=$RESULTS/braker/species --hints=$RESULTS/braker/hintsfile.gff --prot_seq=$DATA/proteins.fa --prg=gth --softmasking --skipAllTraining --cores=12 --AUGUSTUS_CONFIG_PATH=/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/Augustus/config 
