#!/bin/bash

#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=GWAS_GEMMA
#SBATCH --time=0-01:00:00

# Extract all arguments except the last one
args_without_last="${@:1:$#-1}"

# Extract the last argument
last_arg="${!#}"

Method=$last_arg
if [ "$Method" == "minimap2" ]; then
    ./homology.sh $args_without_last
elif [ "$Method" == "blast" ]; then
    ./blast_homology.sh $args_without_last
elif [ "$Method" == "mummer" ]; then
    ./homology_mummer.sh $args_without_last
else
    ./blast_homology.sh $args_without_last; homology.sh $args_without_last 
fi

# USAGE: ./master_script.sh Mechanitis_messenoides Melinaea_menophilus Hypothyris_anastasia (empty (both methods)|minimap2|blast|mummer)
