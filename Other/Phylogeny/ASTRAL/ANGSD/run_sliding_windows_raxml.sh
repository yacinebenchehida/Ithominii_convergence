#!/bin/bash

#SBATCH --mem=20GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-24:00:00
#SBATCH --array=1-21

#####################
# Load used modules #
#####################
module load RAxML/8.2.12-GCC-10.2.0-pthreads-avx2
module load Biopython/1.81-foss-2022b
PXRR="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/phyx/src/pxrr"

######################################
# Select the scaffold to be analysed #
######################################
SCAFFOLD=$(cat /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/ASTRAL/Inputs/list_scaffold.txt| awk -v scaffold=${SLURM_ARRAY_TASK_ID} 'NR==scaffold')

#########################################################################
# Split each scaffold into windows of 5 kb and run raxml on each window #
#########################################################################
./sliding_windows_raxml.py ../Results/"$SCAFFOLD"/"$SCAFFOLD"_combined.fasta ../Results/"$SCAFFOLD"

#####################################################################
# Combine the inferred phylogenies for all windows in a single file #
#####################################################################
cat ../Results/"$SCAFFOLD"/RAxML_bipartitionsBranchLabels* >> ../Results/Combine_phylogeny.txt
cat ../Results/Combine_phylogeny.txt |perl -pe 's/scaffold_(\d+)_//g' >../Results/curated_combined_phylogeny.txt
rm ../Results/Combine_phylogeny.txt
