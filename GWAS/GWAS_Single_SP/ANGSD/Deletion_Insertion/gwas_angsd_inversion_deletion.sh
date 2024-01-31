#! /bin/bash

#SBATCH --job-name=ANGGWAS
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-01:00:00

########## Load modules
module load SAMtools/1.17-GCC-12.2.0
module load BCFtools/1.15.1-GCC-11.3.0
module load HTSlib/1.15.1-GCC-11.3.0

########## Set used paths
ANGSD="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/angsd/angsd"
BAM_FOLDER="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Results/1_sorted_dedup_bam"
WD="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/ANGSD/Inputs/Deletion_Insertion"
REF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/Melinaea_menophilus/CAM015033_Melinaea_menophilus_ssp_nov.fa"

########## Create an empty file with list of the path to the bam files
touch $WD/bam_list.txt

########## Create a bam file for each sample for the region around the gap 
cd $BAM_FOLDER
cat $WD/menophilus_list.txt|while read line
do
	samtools view -b -o $WD/SUPER_18_"$line"_gap.bam $line/"$line"*.bam SUPER_18:14819397-15823591
	readlink -f $WD/*$line*bam >> $WD/bam_list.txt
done

######## Call SNPs using angsd
cd $WD
$ANGSD -b $WD/bam_list.txt -dobcf 1 -gl 2 -dopost 1 -domaf 2 -domajorminor 1 -docounts -doGeno 4 -dosnpstat -snp_pval 0.7 -out $WD/GWAS_around_peak_menophilus

######## Visualise the genotypes with bcftools
bcftools view $WD/GWAS_around_peak_menophilus.bcf|bcftools query -f '%POS  [ %GT]\n'] 

######## Create phenotype file
cat $WD/menophilus_list.txt |while read line; do grep $line /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Inputs/Melinaea_menophilus/GEMMA_encoding_phenotype_Melinaea_menophilus.txt; done |awk '{print $3}'|perl -pe 's/3/1/g'|perl -pe 's/2/0/g'|perl -pe 's/-9/-999/g' > $WD/phenotype_file.txt

######## Run angsd based GWAS
$ANGSD -yQuant $WD/phenotype_file.txt -doAsso 2 -GL 2 -doPost 1 -out GWAS_around_peak_menophilus -doMajorMinor 1 -SNP_pval 0.7 -doMaf 1 -Pvalue 1 -ref $REF  -bam $WD/bam_list.txt 

######## Plot Gwas results
zcat $WD/GWAS_around_peak_menophilus.lrt0.gz |grep -v -E "nan|Chromosome"| awk '$8 < 1'|awk '{print $1"\t"$2"\t"$8}' > $WD/GWAS_around_peak_menophilus.txt
Rscript Plot.R $WD/GWAS_around_peak_menophilus.txt $WD/scaffold_order.txt 136000
