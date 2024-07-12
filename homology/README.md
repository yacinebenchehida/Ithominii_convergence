# GWAS peaks alignments

This folder contains all the scripts used to perform the peaks alignment in sliding windows using nucmer from the MUMmer package.

## 0) Running the whole pipeline

The whole pipeline can be run using the command below: 

``` bash
./master.sh species1 species2 cortex mummer
```
This pipeline will run the analysis and generate an "alignment plot" for the two species of interest around the cortex region.  
It requires  [MUMmer](https://mummer.sourceforge.net/manual/) and [Biopython](http://biopython.org/) to work. 

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
