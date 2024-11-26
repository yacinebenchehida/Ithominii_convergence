#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=transP

module load Pysam/0.22.0-GCC-12.3.0

VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/multisp/Melinaea_menophilus/Multisp_menophilus.DP4_GQ5_QUAL5_invariants.vcf.gz"
POP="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/Transpolymorphism/Data/yellow_band_species.txt"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/Transpolymorphism/Results/SUPER_4_menophilus"

mkdir -p $RESULTS

#python3 transP_v3.py  -v $VCF -p $POP -c SUPER_18 -o SUPER_18_Menophilus -w 50000 -d ../Results --highlight_start 15186353 --highlight_end 15313931
#python3 transP_v3.py  -v $VCF -p $POP -c SUPER_5 -o SUPER_5_Menophilus -w 50000 -d $RESULTS --highlight_start 25892884 --highlight_end 26028482
python3 transP_v3.py  -v $VCF -p $POP -c SUPER_4 -o SUPER_4_Menophilus -w 50000 -d $RESULTS --highlight_start 16235447 --highlight_end 16392008
