#!/bin/bash
#Author: Edd Page

#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=homer

INPUTS="/mnt/scratch/projects/biol-specgen-2018/edd/Transcription/ROIs/Results"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/edd/Transcription/homer/Results/"

rm $INPUTS/no_YB.fasta
rm $INPUTS/YB.fasta

cat $INPUTS/deceptus/* >> $INPUTS/no_YB.fasta
cat $INPUTS/phasianita/* >> $INPUTS/no_YB.fasta

cat $INPUTS/messenoides/* >> $INPUTS/YB.fasta


mkdir -p $RESULTS/mechanitis
rm -r $RESULTS/mechanitis/*

findMotifs.pl $INPUTS/YB.fasta mechanitis $RESULTS/mechanitis/ -mset insects -fastaBg $INPUTS/no_YB.fasta

sed 's/.*>.*//' $RESULTS/mechanitis/homerMotifs.all.motifs > $RESULTS/mechanitis/mechanitis_vertical.txt
