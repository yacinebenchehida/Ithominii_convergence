# Genotype maps

This folder contains all the scripts used to generate the genotype maps (Extended Data Fig. 3-9,11-13,15)

## 0) Prepare the data for the analysis

The first step uses PLINK to convert the raw input file into bed, bim and fam format files. These different files used by GEMMA to perform downstream analyses. 

``` bash
$PLINK --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --pheno $PHENOTYPES --make-bed --out $RESULTS
```

## 1) Get pairwise relatedness between individuals

To control for  population structure, GEMMA computed the relatedness among samples. We used the following command to compute the pairwise relastedness matrix:

``` bash
$GEMMA -bfile $RESULTS_FOLDER -gk 1 > gemma_relatedness2.out
```

## 2) Run the GWAS

GEMMA was ran using the Univariate Linear Mixed Model. We employed the log-likelihood ratio test if the size effect of each SNP was significantly different from 0. 

``` bash
$GEMMA -bfile $RESULTS_FOLDER -lmm 4 -o assoc.gemma -k relatedness.cXX.txt > gemma.out
``` 

## 3) Plot the Manhattan plot

Manhattan plot summarising the results along the whole genome were plotted using R.

- In order to scaffolds plotted by increasing order we first ran the scaffold_size.py script on the reference genome.
``` bash
# Create file with scaffold sorted by ascending size
module load  Biopython/1.79-foss-2022a
python ./scaffold_size.py  $REFERENCE|sort -k 2 -nr|awk '{print $1}' > scaffold_order.txt
```

- The to get the number of SNPs used for the analyses (used to define the bonferroni threshold) we defined the THRESHOLD variable:
``` bash
THRESHOLD=$(cat $RESULTS/*assoc.gemma.log.txt|grep "analyzed SNPs/var" |awk '{print $7}')
``` 

- Finally the GWAS Manhattan plot is plotted using the R script Plot.R.
``` bash
# Make a png plot of the gwas 
Rscript ./Plot.R $RESULTS/*.assoc.gemma.assoc.txt $RESULTS/scaffold_order.txt $THRESHOLD $RESULTS
```


