# RELATE

This folder contains all the scripts used to perform the RELATE analyses.

## 0) Running the whole pipeline

The whole pipeline can be run using the command below: 

``` bash
sbatch ./Relate_launcher \
-v VCF_FILE \
-c CHROMOSOME/SCAFFOLD \
-s START_POSITION \ 
-e END_POSITION \ 
-f POPULATION_FILE \
-t TAXA1,TAXA2,TAXA3,TAXA4 \ # List of taxa to be included in the analysis (separated by a comma)
--snps SNP1,SNP2,SNP3,SNP4 \ # List of SNP positions for which a pdf tree will be generated (separated by a comma)
-r ANCESTRAL_STATE \ # Taxa used to set the ancestral state (This taxa should also be included in the -t flag)
-o OUTPUT_PATH \
-n OUTPUT_NAME
```
This pipeline will run the analysis and generate an "alignment plot" for the two species of interest around the cortex region.  
It requires  [SHAPEIT4](https://odelaneau.github.io/shapeit4/), [VCFTOOLS](https://vcftools.github.io), [BCFTOOLS](https://samtools.github.io/bcftools/) and [RELATE](https://myersgroup.github.io/relate/) to work. 

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







