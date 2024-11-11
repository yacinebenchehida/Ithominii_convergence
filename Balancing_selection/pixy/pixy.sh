#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=pixy

module load BCFtools/1.19-GCC-13.2.0

POP_FILE="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/pixy/Inputs/$2"
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/multisp/Melinaea_mothone/multisp_Melinaea_mothone_genotypeGVCF.intervals_15.filters.DP4_GQ5_QUAL5.invariants.vcf.gz"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/pixy/Results/$1"
CHR="SUPER_3"
START="32942064"
END="39553786"
PEAK_START="39553686"
PEAK_END="39553786"
mkdir -p $RESULTS

awk -v chr="$CHR" -v start="$START" -v end="$END" '$1 == chr && $2 >= start && $3 <= end' /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/pixy/Inputs/masked.txt > $RESULTS/region_to_remove.txt
mask="$RESULTS/region_to_remove.txt"

bcftools view --regions $CHR:$START-$END -S <(cat $POP_FILE|awk '{print $1}') --threads 8 $VCF -Oz -o $RESULTS/"$1".vcf.gz
bcftools view  -T "^$mask" $RESULTS/"$1".vcf.gz -o $RESULTS/"$1"_filtered.vcf.gz
tabix $RESULTS/"$1"_filtered.vcf.gz

NEW_VCF=$RESULTS/"$1"_filtered.vcf.gz
POPULATION_FILE_NO_MISS="$RESULTS/population_file_no_miss_data.txt"
cat $POP_FILE | grep -v "\-9" > $POPULATION_FILE_NO_MISS
NEW_PHENOTYPE=$POPULATION_FILE_NO_MISS

pixy --stats pi \
--vcf $NEW_VCF \
--populations $NEW_PHENOTYPE \
--window_size 10000 \
--n_cores 8 \
--chromosomes $CHR \
--output_folder $RESULTS \
--output_prefix $1

cat $RESULTS/"$1"_pi.txt| grep -v "NA" | awk '{printf "%.0f\t%s\n", ($3+$4)/2, $5}' > $RESULTS/"$1"_pi2plot.txt

module load  R/4.2.1-foss-2022a
Rscript ./plot_nucldiv.R $RESULTS/"$1"_pi2plot.txt $PEAK_START $PEAK_END $RESULTS
mv $RESULTS/nucleotide_diversity.pdf $RESULTS/"$1"_nucleotide_diversity.pdf
