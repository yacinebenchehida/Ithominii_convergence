#! /bin/bash

#SBATCH --mem=25GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=GWAS_GEMMA
#SBATCH --time=0-02:00:00

module load R/4.2.1-foss-2022a
module load Pillow/9.4.0-GCCcore-12.2.0
module load BCFtools/1.15.1-GCC-11.3.0

# Generate plots
Rscript New_Combine_GWAS_plots.R $@  2> warnings.txt

# Crop png
#for i in $@; do echo $i; python ./crop_png.py "$i"_Zoom_GWAS.png Zoom_GWAS_"$i"; done

# Crop pdf

module purge
module load PyPDF2/1.26.0-foss-2019b-Python-2.7.16

# Get all command line arguments into an array
args=("$@")

# Determine the index of the last element
last_index=$((${#args[@]} - 1))

# Loop over all elements except the last one
for ((i = 0; i < last_index; i++)); do
	echo ${args[i]}; python3 ./cropper.py "${args[i]}"_Zoom_GWAS.pdf Zoom_GWAS_"${args[i]}".pdf
done

# Move results to results folder 
#rm  "$i"_Zoom_GWAS.png  #"$i"_Zoom_GWAS_features.pdf 
mv *png *pdf ../Results/"${args[last_index]}"
