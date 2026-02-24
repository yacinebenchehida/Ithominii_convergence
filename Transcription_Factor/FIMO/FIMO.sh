#!/bin/bash
#Author: Edd Page

#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=FIMO

module load MEME/5.5.0-gompi-2022b

INPUTS="/mnt/scratch/projects/biol-specgen-2018/edd/Transcription/FIMO/Inputs"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/edd/Transcription/FIMO/Results"
ROIs="/mnt/scratch/projects/biol-specgen-2018/edd/Transcription/ROIs/Results"

#Mechanitis combined:

fasta-get-markov -m 5 $INPUTS/Mechanitis_SUPER_6.fasta > $INPUTS/Mechanitis_SUPER_6_5th_order_bg_model.txt

fimo -thresh 0.0001 --bfile $INPUTS/background/Mechanitis_SUPER_6_5th_order_bg_model.txt -oc $RESULTS/mechanitis_combined $INPUTS/all_insect_motifs_JASPAR.txt $INPUTS/mechanitis_SUPER_6:6877302-6878798_combined.fasta

fimo -thresh 0.0001 --bfile $INPUTS/background/Mechanitis_SUPER_6_5th_order_bg_model.txt -oc $RESULTS/mechanitis_bicolora $INPUTS/all_insect_motifs_JASPAR.txt $INPUTS/mechanitis_SUPER_6:6877302-6878798_combined.fasta
