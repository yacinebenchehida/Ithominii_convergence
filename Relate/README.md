# RELATE

This folder contains all the scripts used to perform the RELATE analyses.

## 0) Running the whole pipeline

The whole pipeline can be run using the command below: 

``` bash
sbatch ./Relate_launcher -v /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/multisp/Melinaea_marsaeus/multisp_Melinaea_marsaeus_genotypeGVCF.intervals_8.filters.DP4_GQ5_QUAL5.invariants.vcf.gz -c SUPER_2 -s 23556241 -e 27031063 -f ../Inputs/Input_2_use.txt -t mothone,outgroup,isocomma,satevis,marsaeus_rileyi,marsaeus_phasiana,lilis,flavomacula --snps 25973722,25973735,25973747,25976514,25976548,25976654,25993654,25993663,25993709,25996282 -r outgroup -o ../Results/Melinaea_marsaeus -n mothone_ludovica_isocomma_satevis_marsaeus_rileyi_marsaeus_phasiana_lilis_flavomacula
```
This pipeline will run the analysis and generate an "alignment plot" for the two species of interest around the cortex region.  
It requires  [SHAPEIT4](https://odelaneau.github.io/shapeit4/), [VCFTOOLS](https://vcftools.github.io), [BCFTOOLS](https://samtools.github.io/bcftools/) and [RELATE](https://myersgroup.github.io/relate/) to work. 

## 1) Extracting the regions to align from the reference genome

The command below extracts the genomics regions we want to align from the two species we want to compare. This step uses the python script selection_seq_interval.py. 

``` bash
python selection_seq_interval.py  fasta_reference_species1.fasta Scaffold_name starting_position_species1 end_position_species1 Species1_name > Species1_name.fasta
python selection_seq_interval.py  fasta_reference_species2.fasta Scaffold_name starting_position_species2 end_position_species2 Species2_name > Species2_name.fasta

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

# Last incomplete windows (smaller than 1000bp)
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

## 3) Plot the alignment plots

Plots summarising the results along the regions of interest were plotted in R using the script plotting_nucmer_windows.R.

``` bash
Rscript ./plotting_nucmer_windows.R mapping_nucmer_sliding_windows_sp1_sp2.txt sp1 sp2
```
