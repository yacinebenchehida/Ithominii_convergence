# RELATE (WORK IN PROGRESS)

This folder contains all the scripts used to perform the RELATE analyses.

## 0) Running the whole pipeline

The whole pipeline can be run using the command below: 

``` bash
./Relate_launcher \
-v VCF_FILE \
-c CHROMOSOME/SCAFFOLD \
-s START_POSITION \ 
-e END_POSITION \ 
-f POPULATION_FILE \
-t TAXA1,TAXA2,TAXA3,TAXA4 \ # List of taxa to be included in the analysis (separated by a comma)
--snps SNP1,SNP2,SNP3,SNP4 \ # List of SNP positions for which a pdf tree will be generated (separated by a comma)
-r ANCESTRAL_STATE_TAXA1,ANCESTRAL_STATE_TAXA2 \ # Taxa used to set the ancestral state (separated by a comma)
--species NAME1,NAME2,NAME3,NAME4 # List of species acronyms used to set the ancestral state for SNPs not available in the taxa specified using -r
-o OUTPUT_PATH \
-n OUTPUT_NAME

# An example:
/Relate_launcher
-v My_VCF.gz
-c Chr2
-s 24000000
-e 26096070
-f hindwing_black.txt
-t satevis,marsaeus_phasiana,marsaeus_rileyi,mothone
--snps 24000000,25973722,25973735,25973747,25976514,25976548,25976654,25993654,25993663,25993709,25996282,26096070
-r outgroup
--species tarapotensis,satevis,lilis,isocomma,idae,marsaeus,menophilus,mothone,flavo
-o ../Relate_results
-n Example
```
This pipeline will run the analysis and generate an "alignment plot" for the two species of interest around the cortex region.  
It requires  [SHAPEIT4](https://odelaneau.github.io/shapeit4/), [VCFTOOLS](https://vcftools.github.io), [BCFTOOLS](https://samtools.github.io/bcftools/), [RELATE](https://myersgroup.github.io/relate/) [Biopython](https://biopython.org), and [TWISST](https://github.com/simonhmartin/twisst) to work. 

## 1) Population file

The population file should be a tab delimited file containing at least 2 columns: samples names in the VCF and group. If more than 2 columns are provided the penultimate column will be used. The population file looks like this:

``` bash
sample1	group1
sample2 group1
sample3	group2
sample4	group3
...

```

## 2) Recombination rate, mutation rate and effective population size

## 3) Rational behind the RELATE analyses

In this study, the RELATE analyses were performed with four taxa, including an outgroup to define the ancestral allele for each bi-allelic SNP. These analyses, aimed at detecting evidence of introgression at the top GWAS SNPs, were conducted across various taxa combinations following the rationale of the ABBA-BABA test (Durand et al. 2011). The taxa comprised two morphologically distinct sister species (sp1 and sp2), a more distant species (sp3) sharing the wing phenotype of sp1, and an outgroup. Each analysis produced a marginal tree for every SNP that exceeded the significance threshold in the corresponding GWAS.

## 4) Running the pipeline

The pipeline is divided into 8 steps.

- Step 1: Get the phenotype file

``` bash
mkdir -p $OUTPUT_PATH/$PREFIX
if [ -f  $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt ]; then
    rm  $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt
fi
touch  $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt

SPECIES=$(echo $TAXA|perl -pe 's/,/ /g')

for i in $SPECIES; do cat ${POP}| grep $i >>  $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt; done
echo POP FILE READY
```

- Step 2: Get VCF for specific region/individuals
``` bash
bcftools view --regions $CHR:$START-$END -S <(cat $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt|awk '{print $1}') $VCF |  bcftools view -v snps -O v  > $OUTPUT_PATH/$PREFIX/"$PREFIX".vcf
bgzip $OUTPUT_PATH/$PREFIX/"$PREFIX".vcf
tabix $OUTPUT_PATH/$PREFIX/"$PREFIX".vcf.gz
VCF="$OUTPUT_PATH/$PREFIX/$PREFIX.vcf.gz"
echo "SUBSETTED VCF FILE READY"
```

- step 3: phase vcf 

``` bash
echo STARTING PHASING
$SHAPEIT --input $OUTPUT_PATH/$PREFIX/"$PREFIX".vcf.gz --region $CHR --output $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased.vcf 
bgzip $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased.vcf
tabix $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased.vcf.gz
echo PHASING PERFORMED

VCF="$OUTPUT_PATH/$PREFIX/$PREFIX_phased.vcf"
rm $OUTPUT_PATH/$PREFIX/"$PREFIX".vcf.gz
```
- step 4: Convert VCF to relate input files (generate the *haps and *sample)

``` bash
# Create VCF with human chromosome names
echo -e $CHR"\t"1 > $OUTPUT_PATH/$PREFIX/chrom_names.txt
bcftools annotate --rename-chrs $OUTPUT_PATH/$PREFIX/chrom_names.txt $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased.vcf.gz -o $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased_renamed.vcf.gz
rm $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased.vcf.gz
echo VCF RENAMED

# Create a plink linkage map
#vcftools --gzvcf $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased_renamed.vcf.gz --out $OUTPUT_PATH/$PREFIX/"$PREFIX"_linkage_map --plink
#echo PLINK LINKAGE MAP GENERATED

# Create a Relate linkage map from the plink linkage map
#Rscript ./makeRelateMap.r $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased_renamed.vcf.gz $OUTPUT_PATH/$PREFIX/"$PREFIX"_linkage_map.map
Rscript ./createuniformrecmap.r $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased_renamed.vcf.gz LM.bed
echo RELATE LINKAGE MAP GENERATED

# Create inputs required by Relate
$RELATE/bin/RelateFileFormats \
                 --mode ConvertFromVcf \
                 --haps  $OUTPUT_PATH/$PREFIX/"$PREFIX".haps \
                 --sample  $OUTPUT_PATH/$PREFIX/"$PREFIX".sample \
                 -i $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased_renamed
echo HAPS AND SAMPLE FILE GENERATED
```

- step 5: Use outgroup to set ancestral and derived alleles and flip

``` bash
# Find name outgroups and write --indv command
indv_command=$(for i in $(cat $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt|grep $OUTGROUP|awk '{print $1}'); do echo -e "--indv $i"; done |perl -pe 's/\n/ /g')
echo $indv_command

# Get the allele frequencies of the outgroups 
vcftools --gzvcf $OUTPUT_PATH/$PREFIX/"$PREFIX"_phased_renamed.vcf.gz  --freq $indv_command --out $OUTPUT_PATH/$PREFIX/outgroups
echo ALLELE FREQUENCIES ESTIMATED FOR OUTGROUPS

# Get ancestral alleles
awk '{split($5,ref,":"); split($6,alt,":"); 
     if(ref[2]>0.5) print ref[1]; else print alt[1]}' \
     <(grep -v CHROM $OUTPUT_PATH/$PREFIX/outgroups.frq) > $OUTPUT_PATH/$PREFIX/ancestral.alleles

#Swap the ancestral state in the haps file 
awk 'FNR==NR {x2[NR] = $1; next}
  {if(x2[FNR]==$5){ref=$5; alt=$4; $4=ref; $5=alt; 
  for(i=6;i<=NF;i++) if($i==0) $i=1; else $i=0}; print}' \
  $OUTPUT_PATH/$PREFIX/ancestral.alleles $OUTPUT_PATH/$PREFIX/"$PREFIX".haps > $OUTPUT_PATH/$PREFIX/"$PREFIX"_ancestral_state.haps
```

- step 6: Generate the pop label file with the sex column
  
``` bash
(echo -e sample"\t"population"\t"group"\t"sex; awk '{print $0"\t"0}'  $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt) > $OUTPUT_PATH/$PREFIX/"$PREFIX"_relate.poplabels
echo RELATE POPLABELS READY
```

step 7: Run RELATE

``` bash
cd  $OUTPUT_PATH/$PREFIX

$RELATE/bin/Relate \
                 --mode All \
                 -m 2.9e-9 \
                 -N 20000000 \
                 --haps "$PREFIX"_ancestral_state.haps \
                 --sample "$PREFIX".sample \
                 --map "$PREFIX"*plink.map \
                 --seed 384 -o $PREFIX
echo RELATE RUN
```

step 8: Extract results for the SNPs of interest

``` bash
# Get list of SNPs
SNP=$(echo $SNPS|perl -pe 's/,/ /g')
# Initialize min and max with an unset value
min_value=
max_value=

# Iterate through the list of snps positions to find minimum and maximum
for num in $SNP; do
# If min_value is unset or current number is less than min_value, update min_value
    if [[ -z "$min_value" || "$num" -lt "$min_value" ]]; then
        min_value="$num"
    fi

    # If max_value is unset or current number is greater than max_value, update max_value
    if [[ -z "$max_value" || "$num" -gt "$max_value" ]]; then
        max_value="$num"
    fi
done

# Extract Relate results between the minimum and maximum SNPs position
$RELATE/bin/RelateExtract \
                 --mode AncToNewick \
                 --haps "$PREFIX"_ancestral_state.haps \
                 --sample "$PREFIX".sample \
                 --anc  "$PREFIX".anc \
                 --mut "$PREFIX".mut \
                 --first_bp $min_value\
                 --last_bp $max_value\
                 --poplabels "$PREFIX"_relate.poplabels \
                 -o "$PREFIX"_"$CHR"_"$min_value"-"$max_value"

echo RELATE RESULTS EXTRACTED

# Plot Tree for the specified SNPs
for i in $SNP; do \
$RELATE/scripts/TreeView/TreeView.sh \
                --haps "$PREFIX"_ancestral_state.haps \
                --sample "$PREFIX".sample \
                --anc "$PREFIX".anc \
                --mut "$PREFIX".mut \
                --poplabels "$PREFIX"_relate.poplabels \
                --bp_of_interest $i \
                --years_per_gen 1 \
		-o "$PREFIX"_"$i"
done
```
