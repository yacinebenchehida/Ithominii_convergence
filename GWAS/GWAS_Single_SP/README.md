# GWAS analyses

This folder contains all the scripts used to perform the GWAS analyses (GEMMA.sh). The GWAS were performed with the GEMMA software. The script

## 0) Prepare the data for the analysis

The first step (busco.sh) consists in converting the plink input file into a bed file. This step uses PLINK:

``` bash
$PLINK --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --pheno $PHENOTYPES --make-bed --out $RESULTS
```

## 1) Get pairwise relatedness between individuals
``` bash
$GEMMA -bfile $RESULTS_FOLDER -gk 1 > gemma_relatedness2.out
```

## 2) Run the GWAS
``` bash
$GEMMA -bfile $RESULTS_FOLDER -lmm 4 -o assoc.gemma -k relatedness.cXX.txt > gemma.out
``` 

## 3) Plot the Manhattan plot

``` bash
# Create file with scaffold sorted by ascending size
module load  Biopython/1.79-foss-2022a
python $SCAFFOLD_SIZE_SCRIPT  $REF/$FASTA_REF|sort -k 2 -nr|awk '{print $1}' > $RESULTS/scaffold_order_$2.txt

# Get the number of SNPs used for the analyses (used to define the bonferroni threshold)
THRESHOLD=$(cat $RESULTS/*assoc.gemma.log.txt|grep "analyzed SNPs/var" |awk '{print $7}')

# Make a png plot of the gwas (low quality)
Rscript $PLOT_SCRIPT $RESULTS/"$2".assoc.gemma.assoc.txt $RESULTS/scaffold_order_$2.txt $THRESHOLD $RESULTS
```


## 4) Predict genes in the GWAS peaks
``` bash

``` 

