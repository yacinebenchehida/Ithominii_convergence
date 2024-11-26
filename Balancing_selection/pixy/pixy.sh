#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=pixy

# Load modules 
module load BCFtools/1.19-GCC-13.2.0

# Define paths and variables
POP_FILE="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/pixy/Inputs/$2"
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/multisp/Melinaea_menophilus/Multisp_menophilus.DP4_GQ5_QUAL5_invariants.vcf.gz"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/pixy/Results/$1"
CHR="SUPER_4"
START="1"
END="34424040"
PEAK_START="16235447"
PEAK_END="16392008"

# Create working directory
mkdir -p $RESULTS

# Extract regions to mask
awk -v chr="$CHR" -v start="$START" -v end="$END" '$1 == chr && $2 >= start && $3 <= end' /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/pixy/Inputs/masked.txt > $RESULTS/region_to_remove.txt
mask="$RESULTS/region_to_remove.txt"

# Create a pop file with only the individuals we wish to analyse
POPULATION_FILE_NO_MISS="$RESULTS/population_file_no_miss_data.txt"
cat $POP_FILE | grep -v "\-9" > $POPULATION_FILE_NO_MISS
NEW_PHENOTYPE=$POPULATION_FILE_NO_MISS

# Generate a vcf without the masked regions
bcftools view --regions $CHR:$START-$END -S <(cat $NEW_PHENOTYPE|awk '{print $1}') --threads 8 $VCF -Oz -o $RESULTS/"$1".vcf.gz
#tabix $RESULTS/"$1".vcf.gz
bcftools view  -T "^$mask" $RESULTS/"$1".vcf.gz -o $RESULTS/"$1"_filtered.vcf.gz
tabix $RESULTS/"$1"_filtered.vcf.gz
NEW_VCF=$RESULTS/"$1"_filtered.vcf.gz

# Run pixy
pixy --stats pi \
--vcf $NEW_VCF \
--populations $NEW_PHENOTYPE \
--window_size 50000 \
--n_cores 8 \
--chromosomes $CHR \
--output_folder $RESULTS \
--output_prefix $1

# Prepare pixy's output for plotting
cat $RESULTS/"$1"_pi.txt| grep -v "NA" | awk '{printf "%.0f\t%s\n", ($3+$4)/2, $5}' > $RESULTS/"$1"_pi2plot.txt

# Plot results
module load R/4.2.1-foss-2022a
Rscript ./plot_nucldiv.R $RESULTS/"$1"_pi2plot.txt $PEAK_START $PEAK_END $RESULTS
mv $RESULTS/nucleotide_diversity.pdf $RESULTS/"$1"_nucleotide_diversity.pdf
