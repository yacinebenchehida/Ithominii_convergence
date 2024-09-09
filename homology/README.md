# GWAS peaks alignments

This folder contains all the scripts used to perform the peaks alignment in sliding windows using nucmer from the MUMmer package.

## 0) Running the whole pipeline

The whole pipeline can be run using the command below: 

``` bash
./master.sh species1 species2 species3 species4 cortex mummer
```

The plots present in Figure 3 can be generated using this command:
``` bash
# Cortex
./master.sh Mechanitis_messenoides Melinaea_menophilus Melinaea_mothone Hypothyris_anastasia Cortex mummer

# Optix
./master.sh Mechanitis_messenoides Melinaea_menophilus Melinaea_marsaeus Hypothyris_anastasia Optix mummer
```

This pipeline will run the analysis and generate an "alignment plot" for the two species of interest around the cortex region.  
It requires  [MUMmer](https://mummer.sourceforge.net/manual/) and [Biopython](http://biopython.org/) to work. 

The pipeline is divided into few blocks of commands:
1) Extracting the peak regions of each species from the reference genome
2) Run Nucmer in sliding windows of 1kb
3) Get the phenotype and genotype at each SNP in the GWAS (used to estimate the squared Spearman coefficients)
4) Plot the results

## 1) Extracting the regions to align from the reference genome

The command below extracts the genomics regions around the peak of two species of interests.

``` bash
python selection_seq_interval.py  fasta_reference_species1.fasta Scaffold_name starting_position_species1 end_position_species1 Species1_name > Species1_name.fasta
python selection_seq_interval.py  fasta_reference_species2.fasta Scaffold_name starting_position_species2 end_position_species2 Species2_name > Species2_name.fasta
```

The starting and ending position are provided in a script looking like this: 
``` bash
Hypothyris_anastasia    Cortex  4327766   4328081
Hypothyris_anastasia    Optix   9419647   9420916
Melinaea_menophilus     Cortex  15215155  15224528
Melinaea_menophilus     Optix   25912446  25922831
Melinaea_mothone        Cortex  1385004   1398720
Melinaea_marsaeus       Optix   25965242  26021063
Mechanitis_messenoides  Cortex  6876850   6878893
Mechanitis_messenoides  Optix   14398077  14472077
```

## 2) Run nucmer in sliding windows of 1000bp
The block below divides the fasta sequence of sp1 into windows of 1000bp (python script selection_seq_interval_bis.py) and align them against the fasta sequence of sp2 using nucmer. 

``` bash
echo -e "start\tend\tquery\tqueryLen\tqueryStart\tqueryEnd\tSubject\tSubjectLen\tSubjectStart\tSubjectEnd\tIdentity" >  mapping_nucmer_sliding_windows_sp1_sp2.txt
# Define the size of the size genomic region that is going to be plotted
SIZE=$(python3 -W ignore scaffold_size.py Species1_name.fasta|awk '{print $2}')

# Define the window slide and the slide
WINDOWS=1000
SLIDE=1000

# Loop for each pair of species in the define scaffold over windows of the same size
for ((i = 1, j = $WINDOWS; i < $SIZE && j < $SIZE; i = j + 1, j=j+$SLIDE)) 
	do
		python selection_seq_interval_bis.py Species1_name.fasta Scaffold_name $i $j > sp1_"$i"_"$j".fasta # get a fasta for a window of 1000bp
		nucmer --mum -c 20 -b 500 -l 10 --maxgap 500 -p tmp_"$i"_"$j" sp2.fasta sp1_"$i"_"$j".fasta # align the window of 1000bp against sp2
		show-coords -rcl tmp_"$i"_"$j".delta > tmp_nucmer_sp1_sp2_"$i"_"$j".txt # Transform the raw mummer output into an output that shows the coordinate
		(cat tmp_nucmer_sp1_sp2_"$i"_"$j".txt |grep -v "====="|awk 'NR> 4'|perl -pe 's/ +/\t/g' |perl -pe 's/^\t//g'|perl -pe 's/\|\t//g'|awk '{print $12"\t"$8"\t"$1"\t"$2"\t"$13"\t"$9"\t"$3"\t"$4"\t"$7}') > tmp3 # extract the useful information for the mummer results and reshape it for plotting
		awk -v i="$i" -v j="$j" '{print i "\t" j "\t" $0}' tmp3 >> mapping_nucmer_sliding_windows_sp1_sp2.txt # Add window number to the final output
done

# Perform the same on the last incomplete windows (smaller than 1000bp)
j=$SIZE
python selection_seq_interval_bis.py Species1_name.fasta Scaffold_name $i $j > sp1_"$i"_"$j".fasta # get a fasta for a window of 1000bp
nucmer --mum -c 20 -b 500 -l 10 --maxgap 500 -p tmp_"$i"_"$j" sp2.fasta sp1_"$i"_"$j".fasta # align the window of 1000bp against sp2
show-coords -rcl tmp_"$i"_"$j".delta > tmp_nucmer_sp1_sp2_"$i"_"$j".txt # Transform the raw mummer output into an output that shows the coordinate
(cat tmp_nucmer_sp1_sp2_"$i"_"$j".txt |grep -v "====="|awk 'NR> 4'|perl -pe 's/ +/\t/g' |perl -pe 's/^\t//g'|perl -pe 's/\|\t//g'|awk '{print $12"\t"$8"\t"$1"\t"$2"\t"$13"\t"$9"\t"$3"\t"$4"\t"$7}') > tmp3 # extract the useful information for the mummer results and reshape it for plotting
awk -v i="$i" -v j="$j" '{print i "\t" j "\t" $0}' tmp3 >> mapping_nucmer_sliding_windows_sp1_sp2.txt # Add window number to the final output
```

Due to the difficulty of aligning distantly related genomics regions, we used the following flags in nucmer:
nucmer --mum -c 20 -b 500 -l 10 --maxgap 500

**--mum**: matches that are unique in both the reference and query.

**-c 20**: Sets the minimum length of a cluster of matches to 20bp.

**-b 500**: Sets the distance an alignment extension will attempt to extend poor scoring regions before giving up to 500bp.

**-l 10**: Sets the minimum length of a single match to 10bp. (For comparison involving Heliconius melpomene this flag was set to 8).

**-maxgap**: Sets the maximum gap between two adjacent matches in a cluster to 500bp. 


## 3) Get the phenotype and genotype at each SNP in the GWAS (used to estimate the squared Spearman coefficients)


## 4) Plot the alignment plots

Plots summarising the results along the regions of interest were plotted in R using the script plotting_nucmer_windows.R.

``` bash
Rscript ./plotting_nucmer_windows.R mapping_nucmer_sliding_windows_sp1_sp2.txt sp1 sp2
```
