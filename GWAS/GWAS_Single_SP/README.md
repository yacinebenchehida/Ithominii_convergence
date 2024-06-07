# GWAS analyses

This folder contains all the scripts used to perform the GWAS analyses (GEMMA.sh). The GWAS were performed with the GEMMA software. The script

## 0) Prepare the data for the analysis

The first step (busco.sh) consists in converting the plink input file into a bed file. This step uses PLINK:

``` bash
$PLINK --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --pheno $PHENOTYPES --make-bed --out $RESULTS
```

## 1) Get pairwise relatedness between individuals
## 2) Run the GWAS
## 3) Plot the Manhattan plot
## 4) Predict genes in the GWAS peaks

