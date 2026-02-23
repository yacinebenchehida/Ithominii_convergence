#!/bin/bash

##########################
# Load necessary modules #
##########################
module load R/4.2.1-foss-2022a        # Load R module for statistical computing and graphics
module load Biopython/1.81-foss-2022b # Load Biopython module for biological computation
module load MUMmer/3.23-GCCcore-9.3.0 # Load MUMmer module for genome alignment
module load BCFtools/1.19-GCC-13.2.0  # Load BCFtools module for working with VCF files

##########################
# Capture the last argument as gene_name #
##########################
gene_name=${@: -1}  # Get the last argument passed to the script as the gene name
species_list=(${@:1:$(($#-1))})  # Get all arguments except the last one (the species list)

echo "Gene name: $gene_name"          # Print the gene name
echo "Species: ${species_list[@]}"    # Print the list of species

##########################
# Set useful directories #
##########################
Results="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Results"  # Directory for results
Inputs="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Inputs"    # Directory for input files
annotation="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Inputs/${gene_name}_annotation.txt"  # Path to annotation file
GWAS_RESULTS_PATH="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/${gene_name}/GWAS"  # Path to GWAS results
VCF_PATH="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results"  # Path to VCF files
PHENOTYPE_PATH="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/${gene_name}/Phenotypes"

############################################################################################
# Find top SNPs in the GWAS and create a fasta based on a windows of 150kb around the peak #
############################################################################################
for i in ${species_list[@]}  # Loop through each species in the species list
do
    echo $i  # Print the current species
    ref_genome="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$i"  # Path to reference genome for current species
    fasta_sequence=$(ls $ref_genome | grep -E "*.fa$|*.fasta$")  # Find fasta files in the reference genome directory
    
    # Extract scaffold and peak positions for the current species
    Scaffold=$(grep $i $annotation | head -n 1 | awk -F"\t" '{print $2}')
    starting_position=$(python3 ./find_min_max.py $annotation $i | awk '{print $1}')
    ending_position=$(python3 ./find_min_max.py $annotation $i | awk '{print $2}')
    
    echo $Scaffold $starting_position $ending_position  # Print scaffold and peak positions
    
    # Create fasta file with sequences from the defined window around the peak
    python selection_seq_interval.py $ref_genome/$fasta_sequence $Scaffold $starting_position $ending_position $i > $Results/"$gene_name"_"$i".fasta
done

echo "FOUND PEAKS"  # Indicate that peak finding is completed

#######################################################################################
# Get a VCF and the phenotype/genotypes of each species for all SNPs in the GWAS peak #
#######################################################################################
for i in ${species_list[@]}  # Loop through each species
do
    echo $i  # Print the current species
    ref_genome="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$i"  # Path to reference genome for current species
    Scaffold=$(grep -e $i $annotation | head -n 1 | awk -F"\t" '{print $2}')
    starting_peak_position=$(grep -e "$i" peaks | grep -e "$gene_name" | awk -F"\t" '{print $3}')
    ending_peak_position=$(grep -e $i peaks|grep -e $gene_name|awk -F"\t" '{print $4}')
    
    echo -e $Scaffold $starting_peak_position $ending_peak_position  # Print scaffold and peak positions

    # Extract p-values for peaks from GWAS results
    awk -v scaffold="$Scaffold" -v start="$starting_peak_position" -v end="$ending_peak_position" \
        '$1 == scaffold && $3 >= start && $3 <= end {print $3"\t"$4}' "$GWAS_RESULTS_PATH/${i}.txt" > "${i}_peak_pvalues.txt"
    echo "PVALUES IN PEAK EXTRACTED"  # Indicate that p-values extraction is completed

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
    bcftools query -f '%POS  [ %GT]\n'] $VCF > "tmp_${i}.txt"
    grep -f <(awk '{print $1}' "${i}_peak_pvalues.txt") "tmp_${i}.txt" >  "Genotype_${i}.txt"
    rm "tmp_${i}.txt"
    
    # Figure out number of SNPs in the whole GWAS output
    TOTAL=$(cat "$GWAS_RESULTS_PATH/${i}.txt" |wc -l)
 
    #  Create a genotype phenotype input for each SNP and each gwas peak 
    if [ -f "${i}_genotype_phenotype_input.txt" ]; then
        rm "${i}_genotype_phenotype_input.txt"
    fi

    cat "Genotype_${i}.txt"|while read line; do
	SNP=$(echo $line|awk '{print $1}')
	pvalues=$(grep $SNP "${i}_peak_pvalues.txt"|awk '{print $2}')
	paste <(for i in $(seq "$NUMB_SAMPLES"); do echo -e $SNP"\t"$pvalues"\t"$TOTAL; done) <(cat $PHENOTYPE) <(echo $line|awk '{$1=""; print $0}'|perl -pe 's/^ //g'|perl -pe 's/ /\n/g') >>  "${i}_genotype_phenotype_input.txt"
    done

    #rm  "Genotype_${i}.txt" $VCF* "${i}_peak_pvalues.txt"
    echo -e "GENOTYPE PHENOTYPE FILE READY FOR  ${i}"
done


########################################
# Run Nucmer in sliding windows of 1kb #
########################################
values=("${species_list[@]}")  # Prepare a list of species for sliding window comparisons

for ((k = 0; k < ${#values[@]} - 1; k += 1)); do # Iterate through adjacent pairs of species
    current="${values[k]}"
    next="${values[k + 1]}"
    echo -e "START SLIDING WINDOWS $current $next"
    echo -e "start\tend\tsubject\tsubjectLen\tsubjectStart\tsubjectEnd\tquery\tqueryLen\tqueryStart\tqueryEnd\tIdentity" > mapping.txt
    
    # Get size and scaffold information for the current species
    SIZE=$(python3 -W ignore scaffold_size.py $Results/"$gene_name"_"$current".fasta | awk '{print $2}')
    SCAFFOLD=$(python3 -W ignore scaffold_size.py $Results/"$gene_name"_"$current".fasta | awk '{print $1}')
    
    # Define window size and slide
    WINDOWS=1000
    SLIDE=1000
    echo $SIZE $current
    
    # Loop over sliding windows for the given scaffold
    for ((i = 1, j = $WINDOWS; i < $SIZE && j < $SIZE; i = j + 1, j = j + $SLIDE))
    do
        echo -E "$i $j $SIZE"
        python selection_seq_interval_bis.py $Results/"$gene_name"_"$current".fasta $SCAFFOLD $i $j > "$current"_"$i"_"$j".fasta
        nucmer --mum -c 20 -b 500 -l 10 --maxgap 500 -p tmp_"$i"_"$j" $Results/"$gene_name"_"$next".fasta "$current"_"$i"_"$j".fasta
        show-coords -rcl tmp_"$i"_"$j".delta > tmp_nucmer_"$current"_"$next"_"$i"_"$j".txt
        (cat tmp_nucmer_"$current"_"$next"_"$i"_"$j".txt | grep -v "=====" | awk 'NR> 4' | perl -pe 's/ +/\t/g' | perl -pe 's/^\t//g' | perl -pe 's/\|\t//g' | awk '{print $12"\t"$8"\t"$1"\t"$2"\t"$13"\t"$9"\t"$3"\t"$4"\t"$7}') > tmp3
        awk -v i="$i" -v j="$j" '{print i "\t" j "\t" $0}' tmp3 >> mapping.txt
        rm "$current"_"$i"_"$j".fasta tmp*
    done

    # Process the last window
    j=$SIZE
    echo -E "$i $j $SIZE"
    python selection_seq_interval_bis.py $Results/"$gene_name"_"$current".fasta $SCAFFOLD $i $j > "$current"_"$i"_"$j".fasta
    nucmer --mum -c 20 -b 500 -l 10 --maxgap 500 -p tmp_"$i"_"$j" $Results/"$gene_name"_"$next".fasta "$current"_"$i"_"$j".fasta
    show-coords -rcl tmp_"$i"_"$j".delta > tmp_nucmer_"$current"_"$next"_"$i"_"$j".txt
    (cat tmp_nucmer_"$current"_"$next"_"$i"_"$j".txt | grep -v "=====" | awk 'NR> 4' | perl -pe 's/ +/\t/g' | perl -pe 's/^\t//g' | perl -pe 's/\|\t//g' | awk '{print $12"\t"$8"\t"$1"\t"$2"\t"$13"\t"$9"\t"$3"\t"$4"\t"$7}') > tmp3
    awk -v i="$i" -v j="$j" '{print i "\t" j "\t" $0}' tmp3 >> mapping.txt
    rm "$current"_"$i"_"$j".fasta tmp*

    # Combine the results
    mkdir -p $Results/mummer/sliding_windows_mapping
    cat mapping.txt | awk 'NF > 2' > $Results/mummer/sliding_windows_mapping/mapping_nucmer_sliding_windows_"$current"_"$next".txt
    rm mapping.txt
    echo -E "END SLIDING WINDOWS NUCMER $current $next"
done

#######################
# plot nucmer results #
#######################
echo "START PLOTTING"  # Indicate that plotting is starting
Rscript ./plotting_nucmer_windows.R  $annotation  # Run R script for plotting
mv *pdf $Results/mummer/sliding_windows_mapping/  # Move resulting PDFs to the results directory
#rm *_genotype_phenotype_input.txt
echo -E "ALL SLIDING WINDOWS FINISHED"  # Indicate that all sliding window analysis is completed
