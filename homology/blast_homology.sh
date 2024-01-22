#!/bin/bash

##########################
# load necessary modules #
##########################
module load R/4.2.1-foss-2022a
module load Biopython/1.81-foss-2022b
module load minimap2/2.26-GCCcore-12.2.0

##########################
# Set useful directories #
##########################
Results="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Results"
data_gwas="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data"
Inputs="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Inputs"

############################################################################################
# Find top SNPs in the GWAS and create a fasta based on a windows of 150kb around the peak #
############################################################################################
for i in $@
do
	ref_genome="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$i"
	fasta_sequence=$(ls $ref_genome|grep -E "*.fa$|*.fasta$")
	cat $data_gwas/"$i".txt |awk '$4 < 0.000000001'|LC_ALL=C sort -g -k 4|head -n 10|awk '{print $1"\t"$3"\t"$4}' > top_peak_positions_"$i".txt	
	values=$(Rscript position_extra.R top_peak_positions_"$i".txt)
	echo $values
	python selection_seq_interval.py  $ref_genome/$fasta_sequence $values $i> $Results/Cortex_"$i".fasta
	#rm top_peak_positions_"$i".txt
done

echo FOUND PEAKS

###########################
# Homology based on blast #
###########################
module purge
module load BLAST+/2.14.0-gompi-2022b
module load Biopython/1.81-foss-2022b

echo START BLAST
values=("$@")

for ((sp = 0; sp < ${#values[@]} -1; sp += 1)) # Iterate through adjacent species pairs
do	
	# Define the pair of species that are going to be analyses in this loop
	current="${values[sp]}"
	next="${values[sp + 1]}"
	echo -e "START SLIDING WINDOWS $current $next"
	makeblastdb -in $Results/Cortex_"$current".fasta -dbtype nucl -input_type fasta -out $current -title $current
	makeblastdb -in $Results/Cortex_"$next".fasta -dbtype nucl -input_type fasta -out $next -title $next
	echo -e "start\tend\tquery\tsubject\tidentity\talignmentlength\tmismatches\tgapopens\tquerystart\tqueryend\tsubjectstart\tsubjectend\tevalue\tbitscore" > mapping_"$current"_"$next".txt

	# Define the size of the size genomic interval and the scaffold name that is going to be plotted
	SIZE=$(python3 -W ignore scaffold_size.py $Results/Cortex_"$next".fasta|awk '{print $2}')
	SCAFFOLD=$(python3 -W ignore scaffold_size.py $Results/Cortex_"$next".fasta|awk '{print $1}')

	# Define the window slide and the slide (here without overlap so both have the same size)
	WINDOWS=1000
	SLIDE=1000
	
	# Loop for each pair of species in the define scaffold over windows of the same size
	for ((i = 1, j = $WINDOWS; i < $SIZE && j < $SIZE; i = j + 1, j=j+$SLIDE)) 
	do
		python3 selection_seq_interval_bis.py $Results/Cortex_"$next".fasta $SCAFFOLD $i $j > "$next"_"$i"_"$j".fasta
		RES=$(blastn -task blastn -query "$next"_"$i"_"$j".fasta -db $current -outfmt 7|grep -v "#"|head -n 1) 
		echo -e $i"\t"$j"\t"$RES |perl -pe 's/ /\t/g' >> mapping_"$current"_"$next".txt
		rm "$next"_"$i"_"$j".fasta
	done
	
	#tblastn cortex in both genomes
	echo -e "start\tend\tquery\tsubject\tidentity\talignmentlength\tmismatches\tgapopens\tquerystart\tqueryend\tsubjectstart\tsubjectend\tevalue\tbitscore" > cortex_blasting.txt
	ref_genome1="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$current"
        fasta_sequence1=$(ls $ref_genome|grep -E "*.fa$|*.fasta$")
	makeblastdb -in $ref_genome1/$fasta_sequence1 -dbtype nucl -input_type fasta -out $current -title $current	
	tblastn -query ../Inputs/cortex.fasta -db "$current" -outfmt 7|grep -v "#"|head -n 1 >> cortex_blasting.txt
	
	ref_genome2="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$next"
        fasta_sequence2=$(ls $ref_genome|grep -E "*.fa$|*.fasta$")
        makeblastdb -in $ref_genome2/$fasta_sequence2 -dbtype nucl -input_type fasta -out $next -title $next
	tblastn -query ../Inputs/cortex.fasta -db "$next" -outfmt 7|grep -v "#"|head -n 1 >> cortex_blasting.txt

	# Combine the results
done
