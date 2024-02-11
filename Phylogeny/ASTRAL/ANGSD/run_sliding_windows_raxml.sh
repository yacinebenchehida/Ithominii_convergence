#!/bin/bash

#SBATCH --mem=8GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-10:00:00
#SBATCH --array=1-25

module load RAxML/8.2.12-GCC-10.2.0-pthreads-avx2
SCAFFOLD=$(cat /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/ASTRAL/Inputs/list_scaffold.txt| awk -v scaffold=${SLURM_ARRAY_TASK_ID} 'NR==scaffold')

./sliding_windows_raxml.py ../Results/"$SCAFFOLD"/"$SCAFFOLD"_combined.fasta ../Results/"$SCAFFOLD"
