#!/bin/bash

# Load bcftools
module load BCFtools/1.15.1-GCC-11.3.0

# Set useful variables
VCF=$1
SCAFFOLD=$2
START=$3
END=$4
GENE=$5
SPECIES=$6
PHENOTYPE="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/${GENE}/Phenotypes/${SPECIES}.txt"
NUMB_SAMPLES=$(cat $PHENOTYPE|wc -l)

# Get a vcf for the peak regions
bcftools view --regions $SCAFFOLD:$START-$END $VCF > "$SPECIES"_portorico

# Get genotype for the peak regions
bcftools query -f '%POS  [ %GT]\n'] "$SPECIES"_portorico > "$SPECIES"_jamaisdeuxsanstrois
rm "$SPECIES"_portorico

# Extract SNPs above the threshold of significance
cat "$SPECIES"_tmp.txt|awk 'NR > 1' > "$SPECIES"_tmp2.txt
rm  "$SPECIES"_tmp.txt
grep -f "$SPECIES"_tmp2.txt "$SPECIES"_jamaisdeuxsanstrois > "$SPECIES"_genotype_top_SNPS.txt
rm "$SPECIES"_tmp2.txt "$SPECIES"_jamaisdeuxsanstrois 

# Make genotype-phenotype input file for R
> "$SPECIES"_genotype_phenotype_input.txt

cat "$SPECIES"_genotype_top_SNPS.txt|while read line; do
	SNP=$(echo $line|awk '{print $1}')
	paste <(for i in $(seq "$NUMB_SAMPLES"); do echo $SNP; done) <(cat $PHENOTYPE) <(echo $line|awk '{$1=""; print $0}'|perl -pe 's/^ //g'|perl -pe 's/ /\n/g') >> "$SPECIES"_genotype_phenotype_input.txt
done

rm "$SPECIES"_genotype_top_SNPS.txt
