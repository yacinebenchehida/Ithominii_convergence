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

All these four steps are implemented in the script homology_mummer.sh.   

## 1) Extracting the regions to align from the reference genome

The command below extracts the genomics regions around the peak of two species of interests.

``` bash
python selection_seq_interval.py  fasta_reference_species1.fasta Scaffold_name starting_position_species1 end_position_species1 Species1_name > Species1_name.fasta
python selection_seq_interval.py  fasta_reference_species2.fasta Scaffold_name starting_position_species2 end_position_species2 Species2_name > Species2_name.fasta
```

The starting and ending position are provided in a script looking like this: 
``` bash
sp1  Gene1  start  end
sp1  Gene1  start  end
sp2  Gene2  start  end
sp2  Gene2  start  end
sp3  Gene3  start  end
sp3  Gene3  start  end
sp4  Gene4  start  end
sp4  Gene4  start  end
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

The  block of commands below was applied to each species. It:
1) Finds the scaffold, start and ending positions for each gwas peak region
2) Generates an indexed VCF file for each gwas peak region
3) Extracts the phenotype of each sample in the VCF
4) Extracts the genotype all SNPs present in the gwas peak region
5) Puts together all the information in a single text file used in R to generate the plot  

``` bash
for i in ${species_list[@]}  # Loop through each species
do
    echo $i  # Print the current species
    ref_genome="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$i"  # Path to reference genome for current species
    Scaffold=$(grep -e $i $annotation | head -n 1 | awk -F"\t" '{print $2}')
    starting_peak_position=$(grep -e "$i" peaks | grep -e "$gene_name" | awk -F"\t" '{print $3}')
    ending_peak_position=$(grep -e $i peaks|grep -e $gene_name|awk -F"\t" '{print $4}')
    
    echo -e $Scaffold $starting_peak_position $ending_peak_position  # Print scaffold and peak positions

    # Check if the VCF file already exists and remove it if it does
    if [ -f "${i}_${Scaffold}_${starting_peak_position}-${ending_peak_position}.vcf.gz" ]; then
        rm "${i}_${Scaffold}_${starting_peak_position}-${ending_peak_position}.vcf.gz"
    fi

    # Extract VCF data for peaks and compress with bgzip and tabix
    bcftools view --regions $Scaffold:$starting_peak_position-$ending_peak_position $VCF_PATH/$i/*.vcf.gz > "${i}_${Scaffold}_${starting_peak_position}-${ending_peak_position}.vcf"
    bgzip "${i}_${Scaffold}_${starting_peak_position}-${ending_peak_position}.vcf"
    tabix "${i}_${Scaffold}_${starting_peak_position}-${ending_peak_position}.vcf.gz"
    VCF="${i}_${Scaffold}_${starting_peak_position}-${ending_peak_position}.vcf.gz"
    echo -e "VCF ready for ${i}"

    # Get genotypes for each samples at each SNP
    PHENOTYPE="${PHENOTYPE_PATH}/${i}.txt"
    NUMB_SAMPLES=$(cat $PHENOTYPE|wc -l)
    bcftools query -f '%POS  [ %GT]\n'] $VCF > "Genotype_${i}.txt"

    #  Create a genotype phenotype input for each SNP and each gwas peak 
    if [ -f "${i}_genotype_phenotype_input.txt" ]; then
        rm "${i}_genotype_phenotype_input.txt"
    fi

    cat "Genotype_${i}.txt"|while read line; do
        SNP=$(echo $line|awk '{print $1}')
        paste <(for i in $(seq "$NUMB_SAMPLES"); do echo $SNP; done) <(cat $PHENOTYPE) <(echo $line|awk '{$1=""; print $0}'|perl -pe 's/^ //g'|perl -pe 's/ /\n/g') >>  "${i}_genotype_phenotype_input.txt"
    done

    #  Remove temporary files
    rm  "Genotype_${i}.txt" $VCF* "${i}_peak_pvalues.txt"
    echo -e "GENOTYPE PHENOTYPE FILE READY FOR  ${i}"
done
```

## 4) Plot the alignment plots

The last setp is to nucmer alignments results along the regions of interest. This step uses R script plotting_nucmer_windows_new.R. It runs like this:

``` bash
Rscript ./plotting_nucmer_windows_new.R annotation_file
```

The annotation file provided looks like this:
``` bash
Mechanitis_messenoides  SUPER_5      Hyd.   like      14389413  14390685
Mechanitis_messenoides  SUPER_5      peak   14398077  14472077  
Mechanitis_messenoides  SUPER_5      Optix  14490502  14494053  
Mechanitis_messenoides  SUPER_5      LRR1   14520006  14526959  
Hypothyris_anastasia    scaffold_17  Hyd.   like      9462002   9463979
Hypothyris_anastasia    scaffold_17  peak   9419647   9420916   
Hypothyris_anastasia    scaffold_17  Optix  9353523   9356174   
Hypothyris_anastasia    scaffold_17  LRR1   9311114   9317579   
Melinaea_menophilus     SUPER_5      Hyd.   like      25896084  25897981
Melinaea_menophilus     SUPER_5      peak   25912446  25922831  
Melinaea_menophilus     SUPER_5      Optix  25994581  25998030  
Melinaea_menophilus     SUPER_5      LRR1   26020125  26026682  
Melinaea_marsaeus       SUPER_2      Hyd.   like      25967145  25968442
Melinaea_marsaeus       SUPER_2      peak   25969242  26021063  
Melinaea_marsaeus       SUPER_2      Optix  26062063  26065512  
Melinaea_marsaeus       SUPER_2      LRR1   26089756  26096070 
```

