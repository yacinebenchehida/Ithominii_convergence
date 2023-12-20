#! /bin/bash

#SBATCH --mem=25GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=GWAS_GEMMA
#SBATCH --time=0-01:00:00

module load R/4.2.1-foss-2022a
module load Pillow/9.4.0-GCCcore-12.2.0

Rscript Combine_GWAS_plots.R $@

for i in $@; do echo $i; python ./crop_png.py "$i"_Zoom_GWAS.png Zoom_GWAS_"$i"; done
rm  "$i"_Zoom_GWAS.png
mv *png ../Results
