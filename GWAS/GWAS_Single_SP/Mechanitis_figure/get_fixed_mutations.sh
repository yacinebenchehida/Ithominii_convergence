#!/bin/bash

# Load the bcftools module, specifying the version (1.15.1 with GCC 11.3.0) for compatibility
module load BCFtools/1.15.1-GCC-11.3.0

# Set variables passed as positional arguments to the script
VCF=$1              # Input VCF file (variant data)
SCAFFOLD=$2         # Chromosome or scaffold identifier
START=$3            # Start position of the region of interest
END=$4              # End position of the region of interest
PHENO=$5            # Phenotype identifier (used to name output files)

# Define the path to the phenotype file for the given phenotype
PHENOTYPE="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_mechanitis/Data/${PHENO}_phenotypes.txt"

# Count the number of samples (lines) in the phenotype file
NUMB_SAMPLES=$(cat $PHENOTYPE | wc -l)

#################################
# Extract VCF data for the region#
#################################

# Filter the VCF file for the specified scaffold and position range, outputting it to a temporary file
bcftools view --regions $SCAFFOLD:$START-$END $VCF > "${PHENO}_portorico"

#########################################
# Extract genotype data for the region #
#########################################

# Use bcftools to output genotypes (GT) for each position in the specified region
# Save the output to another temporary file
bcftools query -f '%POS  [ %GT]\n' "${PHENO}_portorico" > "${PHENO}_jamaisdeuxsanstrois"
rm "${PHENO}_portorico"

################################################
# Filter SNPs based on significance threshold  #
################################################

# Remove the first line from the SNP list file (assuming header) and save it to another temporary file
cat Mechanitis_messenoides_tmp.txt | awk 'NR > 1' > "${PHENO}_tmp2.txt"
rm Mechanitis_messenoides_tmp.txt  # (Optional) Uncomment to delete the original SNP list

# Extract genotypes of SNPs above the threshold by matching SNP positions in the region genotype file
# Save matched genotypes to a file for top SNPs
grep -f "${PHENO}_tmp2.txt" "${PHENO}_jamaisdeuxsanstrois" > "${PHENO}_genotype_top_SNPS.txt"
rm "${PHENO}_tmp2.txt" "${PHENO}_jamaisdeuxsanstrois" 

##############################################
# Create genotype-phenotype input for R      #
##############################################

# Initialize an empty file to store genotype-phenotype information for R input
> "${PHENO}_genotype_phenotype_input.txt"

# For each SNP in the top SNPs file:
cat "${PHENO}_genotype_top_SNPS.txt" | while read line; do
    # Extract the SNP position
    SNP=$(echo $line | awk '{print $1}')
    
    # Write genotype-phenotype data:
    # 1. SNP position repeated for each sample
    # 2. Phenotype values for each sample
    # 3. Genotype information for each sample (split into individual sample genotypes)
    paste <(for i in $(seq "$NUMB_SAMPLES"); do echo $SNP; done) \
          <(cat $PHENOTYPE) \
          <(echo $line | awk '{$1=""; print $0}' | perl -pe 's/^ //g' | perl -pe 's/ /\n/g') \
          >> "${PHENO}_genotype_phenotype_input.txt"
done

# Clean up by removing the top SNPs genotype file
rm "${PHENO}_genotype_top_SNPS.txt"
