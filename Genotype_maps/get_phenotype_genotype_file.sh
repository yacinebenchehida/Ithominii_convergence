#!/bin/bash

##########################
# Load necessary modules #
##########################
module load R
module load BCFtools

########################
# Get script arguments #
########################
SPECIES=$1 
GENE=$2
GWAS=$3
VCF=$4
VCF_multisp=$5
SCAFFOLD=$6
PEAK_START=$7
PEAK_END=$8
PHENOTYPE=$9
NUMB_SAMPLES=$(cat $PHENOTYPE|wc -l)
PHENOTYPE_MULTI=${10}
NUMB_SAMPLES_MULTI=$(cat $PHENOTYPE_MULTI|wc -l)

#echo "Gene name: $GENE"     # Print the gene name
#echo "Species: $SPECIES"    # Print species
#echo "GWAS file: $GWAS"
#echo "VCF: $VCF"
#echo "VCF multispecies: $VCF_multisp"
#echo "Peak position: $SCAFFOLD $PEAK_START $PEAK_END"
echo "Phenotype focal file: $PHENOTYPE contains $NUMB_SAMPLES samples"
echo "Phenotype multi file: $PHENOTYPE_MULTI contains $NUMB_SAMPLES_MULTI samples"

###############################################
# Get SNPs in the GWAS peak for focal species #
###############################################
awk -v scaffold="$SCAFFOLD" -v start="$PEAK_START" -v end="$PEAK_END" \
        '$1 == scaffold && $3 >= start && $3 <= end {print $3"\t"$4}' $GWAS > "${SPECIES}_peak_pvalues.txt"
cat "${SPECIES}_peak_pvalues.txt"|wc -l
echo "PVALUES IN PEAK EXTRACTED"  # Indicate that p-values extraction is completed

###########################################  
# Get a VCF for focal speciees in GWAS peak #
###########################################  
bcftools view --regions $SCAFFOLD:$PEAK_START-$PEAK_END $VCF > "${SPECIES}_focal_peak.vcf"
bgzip "${SPECIES}_focal_peak.vcf"
tabix "${SPECIES}_focal_peak.vcf.gz"

############################################################
# Get genotype for the peak regions from the focal species #
#############################################################
bcftools query -f '%POS  [ %GT]\n'] "${SPECIES}_focal_peak.vcf.gz" > "${SPECIES}_focal_genotypes_tmp.txt"
grep -f <(awk '{print $1}' "${SPECIES}_peak_pvalues.txt") "${SPECIES}_focal_genotypes_tmp.txt" > "${SPECIES}_focal_genotypes.txt"

rm "${SPECIES}_focal_peak.vcf.gz"
rm "${SPECIES}_focal_peak.vcf.gz.tbi"

#####################################################################################################################################
# Function to check if all sample in the VCF are present in the phenotype file + reorder samples in phenotype file in the VCF order #
#####################################################################################################################################
reorder_sample_vcf_order() {
  local VCF=$1
  local PHENOTYPE=$2
  name=$(basename "$PHENOTYPE")
  > reordered_"$name"
  
  # Read sample names from the VCF and check against the phenotype file
 while read -r line; do
    if grep -Pq "$line(\s)+" "$PHENOTYPE"; then
      grep -P "$line(\s)+" "$PHENOTYPE" >> reordered_"$name"
    else
      echo -e "PROBLEM: The sample $line is MISSING from the $PHENOTOYPE file"
      echo -e "Add it to the phenotype file or consider removing it from the VCF"
      echo -e "Quitting"
      echo -e "The problem is the missing sample in the $PHENOTYPE file. Any messge after this line should be ignored"
      exit 1
    fi
  done < <(bcftools query -l "$VCF")
}

########################################################################
# Reorder samples in phenotype file in the same order as the focal VCF #
########################################################################
reorder_sample_vcf_order $VCF $PHENOTYPE
tmp_name=$(basename "$PHENOTYPE")
phenotype=reordered_"$tmp_name"

############################################
# Make genotype-phenotype input file for R #
############################################
> "${SPECIES}_focal_genotype_phenotype_input.txt"

counter=1
cat "${SPECIES}_focal_genotypes.txt" | while read line; do
    SNP=$(echo $line | awk '{print $1}')
    PVALUE=$(sed -n "${counter}p" "${SPECIES}_peak_pvalues.txt"|awk '{print $2}')
    paste <(for i in $(seq "$NUMB_SAMPLES"); do echo $SNP; done) \
          <(cat $phenotype) \
          <(echo $line | awk '{$1=""; print $0}' | perl -pe 's/^ //g' | perl -pe 's/ /\n/g') \
          <(for i in $(seq "$NUMB_SAMPLES"); do echo $PVALUE; done) \
    >> "${SPECIES}_focal_genotype_phenotype_input.txt"
    ((counter++))
done

rm "${SPECIES}_focal_genotypes.txt"
rm reordered_"$tmp_name"
echo "FOCAL SPECIES R INPUT READY"


###########################################  
# Get a VCF for multi species in GWAS peak #
###########################################  
bcftools view --regions $SCAFFOLD:$PEAK_START-$PEAK_END $VCF_multisp > "${SPECIES}_multi_peak.vcf"
bgzip "${SPECIES}_multi_peak.vcf"
tabix "${SPECIES}_multi_peak.vcf.gz"

######################################################
# Get genotype for the peak regions from the multiSP #
######################################################
bcftools query -f '%POS  [ %GT]\n'] "${SPECIES}_multi_peak.vcf.gz" > "${SPECIES}_multi_genotypes.txt"
rm ${SPECIES}_multi_peak.vcf.*

#########################################################################
# Extract SNPs present in the focal species GWAS from the genotype file #
#########################################################################
grep -f <(awk '{print $1}' "${SPECIES}_peak_pvalues.txt") "${SPECIES}_multi_genotypes.txt" > "${SPECIES}_GWAS_SNPS_multisp_genotype.txt"
rm "${SPECIES}_multi_genotypes.txt"

############################################################################################
# Reorder samples in multispecies phenotype file in the same order as the multispecies VCF #
############################################################################################
reorder_sample_vcf_order $VCF_multisp $PHENOTYPE_MULTI
tmp_name=$(basename "$PHENOTYPE_MULTI")
phenotype=reordered_"$tmp_name"

############################################
# Make genotype-phenotype input file for R #
############################################
> "${SPECIES}_multisp_genotype_phenotype_input.txt"

counter=1
cat "${SPECIES}_GWAS_SNPS_multisp_genotype.txt" | while read line; do
    SNP=$(echo $line | awk '{print $1}')
    PVALUE=$(sed -n "${counter}p" "${SPECIES}_peak_pvalues.txt"|awk '{print $2}')
    paste <(for i in $(seq "$NUMB_SAMPLES_MULTI"); do echo $SNP; done) \
          <(cat $phenotype) \
          <(echo $line | awk '{$1=""; print $0}' | perl -pe 's/^ //g' | perl -pe 's/ /\n/g') \
          <(for i in $(seq "$NUMB_SAMPLES_MULTI"); do echo $PVALUE; done) \
    >> "${SPECIES}_multisp_genotype_phenotype_input.txt"
    ((counter++))
done

rm "${SPECIES}_GWAS_SNPS_multisp_genotype.txt" "${SPECIES}_peak_pvalues.txt" reordered_"$tmp_name"
echo "MULTI SPECIES R INPUT READY"
