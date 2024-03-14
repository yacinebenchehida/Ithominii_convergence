#! /bin/bash

#SBATCH --mem=15GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=GWAS_GEMMA
#SBATCH --time=0-02:00:00

module load R/4.2.1-foss-2022a
module load Pillow/9.4.0-GCCcore-12.2.0

# Generate plots
Rscript New_Combine_GWAS_plots.R $@

# Crop png
for i in $@; do echo $i; python ./crop_png.py "$i"_Zoom_GWAS.png Zoom_GWAS_"$i"; done

# Crop pdf

module purge
module load PyPDF2/1.26.0-foss-2019b-Python-2.7.16

for i in $@; do echo $i; python ./crop_pdf.py "$i"_Zoom_GWAS_features.pdf Zoom_GWAS_features_"$i"; done

# Move results to results folder 
rm  "$i"_Zoom_GWAS.png "$i"_Zoom_GWAS_features.pdf 
mv *png *pdf ../Results
