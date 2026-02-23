#! /bin/bash

#SBATCH --mem=40GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=GWAS_GEMMA
#SBATCH --time=0-03:00:00

module load R/4.2.1-foss-2022a

####################
# Set useful paths #
####################
CHEMIN="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Scripts"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Results/$1/$2"
PLINK="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/plink_linux_x86_64_20231018/plink"
GEMMA="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/gemma/gemma-0.98.5-linux-static-AMD64"
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/$1/*.vcf.gz"
REF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$1"
FASTA_REF=$(ls $REF|grep -E "*.fa$|*.fasta$")
PHENOTYPES="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Inputs/$1/GEMMA_encoding_phenotype_"$2".txt"
PLOT_SCRIPT="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Scripts/Plot.R"
SCAFFOLD_SIZE_SCRIPT="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Scripts/scaffold_size.py"
PLOT_HEATMAP="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Scripts/plotting_outlier_heatmap.sh"
PLOT_HEATMAP_R="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Scripts/plotting_outlier_heatmap.R"
PREFIX="$3"

mkdir -p $RESULTS

# Generate bed files with plink 
$PLINK --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --pheno $PHENOTYPES --make-bed --pca --out $RESULTS/$2 --threads 8
cat $RESULTS/"$2".fam|perl -pe 's/2$/0/g'|perl -pe 's/2.5$/0.5/g'|perl -pe 's/3$/1/g' > $RESULTS/tmp.fam
mv $RESULTS/tmp.fam $RESULTS/"$2".fam

# Get pairwise relatedness between individuals
$GEMMA -bfile $RESULTS/$2 -gk 1 -o "$2".relatedness -miss 0.25 -maf 0.1  -outdir $RESULTS > $RESULTS/gemma_"$2"_relatedness2.out

# GWAS 
$GEMMA -bfile $RESULTS/$2 -lmm 4 -o "$2".assoc.gemma -miss 0.25 -maf 0.1  -outdir $RESULTS -k $RESULTS/"$2".relatedness.cXX.txt > $RESULTS/"$2".gemma.out

# Create output for R Manhattan plots
cat  $RESULTS/*.assoc.gemma.assoc.txt| awk '$14 < 0.05 {print $1"\t"$2"\t"$3"\t"$14}'|grep -v "##" > $RESULTS/plot_"$2".txt

# Create file with scaffold sorted by ascending size
module load  Biopython/1.79-foss-2022a
python $SCAFFOLD_SIZE_SCRIPT  $REF/$FASTA_REF|sort -k 2 -nr|awk '{print $1}' > $RESULTS/scaffold_order_$2.txt

# Get the number of SNPs used for the analyses (used to define the bonferroni threshold)
THRESHOLD=$(cat $RESULTS/*assoc.gemma.log.txt|grep "analyzed SNPs/var" |awk '{print $7}')

# Make a png plot of the gwas (low quality)
Rscript $PLOT_SCRIPT $RESULTS/"$2".assoc.gemma.assoc.txt $RESULTS/scaffold_order_$2.txt $THRESHOLD
mv GWAS.png $RESULTS
mv $RESULTS/GWAS.png  $RESULTS/GWAS_"$2".png

# plot SNPs highly associated to difference in phenotype as heatmap for visualization
$PLOT_HEATMAP 	$1 $2 $THRESHOLD > $RESULTS/SNPS_outliers_input.txt
Rscript $PLOT_HEATMAP_R $RESULTS/SNPS_outliers_input.txt /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Inputs/$1/color_code.txt $PREFIX
mv *pdf $RESULTS
