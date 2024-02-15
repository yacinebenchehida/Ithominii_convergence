#!/bin/bash

#SBATCH --mem=5GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=ASTRAL
#SBATCH --time=0-10:50:00

##########################
# load necessary modules #
##########################
module load Java/13.0.2

##############
# Run ASTRAL #
##############
java -Djava.library.path=/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/ASTRAL/lib/ -jar  /mnt/scratch/projects/biol-specgen-2018/yacine/Tools/ASTRAL/astral.5.7.3.jar -i ../Results/curated_combined_phylogeny.txt -o ../Results/astral_results -t 2 2> ../Results/astral_results.logy

