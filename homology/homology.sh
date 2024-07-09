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
data_gwas="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/Cortex/GWAS"
Inputs="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Inputs"
annotation="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Inputs/annotation.txt"

############################################################################################
# Find top SNPs in the GWAS and create a fasta based on a windows of 150kb around the peak #
############################################################################################
for i in $@
do
	ref_genome="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$i"
	fasta_sequence=$(ls $ref_genome|grep -E "*.fa$|*.fasta$")
	#cat $data_gwas/"$i".txt |awk '$4 < 0.000000001'|LC_ALL=C sort -g -k 4|head -n 3|awk '{print $1"\t"$3"\t"$4}' > top_peak_positions_"$i".txt	
	#values=$(Rscript position_extra.R top_peak_positions_"$i".txt)
	Scaffold=$(cat /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/Cortex/annotation_info/annotation.txt|grep $i|head -n 1|awk -F"\t" '{print $2}')
	starting_position=$(python3 ./find_min_max.py $annotation $i|awk '{print $1}')
	ending_position=$(python3 ./find_min_max.py $annotation $i|awk '{print $2}')
	echo $Scaffold $starting_position $ending_position
	python selection_seq_interval.py  $ref_genome/$fasta_sequence $Scaffold $starting_position $ending_position $i> $Results/Cortex_"$i".fasta
	#rm top_peak_positions_"$i".txt
done

echo FOUND PEAKS
#module purge
#module load BLAST+/2.14.0-gompi-2022b
#makeblastdb -in Cortex_"$i".fasta -dbtype nucl -input_type fasta -out $i -title $i

##############################
# Homology based on minimap2 #
##############################
######## 1)based on the asm20 algorithm that compares genome on the bases of a divergence up to 20%
echo START ASM
values=("$@")

for ((i = 0; i < ${#values[@]} -1; i += 1)); do # Iterate through adjacent pairs
	current="${values[i]}"
	next="${values[i + 1]}"
	echo -e "$current\t$next"
	mkdir -p $Results/minimap2/asm20
	echo -e "start\tend\tquery\tqueryLen\tqueryStart\tqueryEnd\tDirection\tSubject\tSubjectLen\tSubjectStart\tSubjectEnd\tMatchingBases\tMatchLen\tQuality" >  $Results/minimap2/asm20/mapping_minimap2_asm20_"$current"_"$next".txt
	minimap2 -x asm20  $Results/Cortex_"$current".fasta  $Results/Cortex_"$next".fasta|cut -f 1-12 >  $Results/minimap2/asm20/mapping_minimap2_asm20_"$current"_"$next".txt
	Rscript plot_syntheny_asm.R $Results/minimap2/asm20/mapping_minimap2_asm20_"$current"_"$next".txt $current $next
	mv *pdf $Results/minimap2/asm20
	echo -e "RAN MINIMAP2 BASED ON ASM20 $current $next"
done
echo ASM DONE

######## 2) based on the reads mapping algorithm of minimap2. This part works in a sliding windows of 1 kb. 
values=("$@")

for ((k = 0; k < ${#values[@]} - 1; k += 1)); do # Iterate through adjacent pairsS
	# Define the pair of species that are going to be analyses in this loop
	current="${values[k]}"
        next="${values[k + 1]}"
	echo -e "START SLIDING WINDOWS $current $next"
	echo -e "start\tend\tquery\tqueryLen\tqueryStart\tqueryEnd\tDirection\tSubject\tSubjectLen\tSubjectStart\tSubjectEnd\tMatchingBases\tMatchLen\tQuality" >  mapping.txt
	# Define the size of the size genomic interval and the scaffold name that is going to be plotted
	SIZE=$(python3 -W ignore scaffold_size.py $Results/Cortex_"$current".fasta|awk '{print $2}')
	SCAFFOLD=$(python3 -W ignore scaffold_size.py $Results/Cortex_"$current".fasta|awk '{print $1}')
	# Define the window slide and the slide (here without overlap so both have the same size)
	WINDOWS=1000
	SLIDE=1000
	echo $SIZE

 	# Loop for each pair of species in the define scaffold over windows of the same size
	for ((i = 1, j = $WINDOWS; i < $SIZE && j < $SIZE; i = j + 1, j=j+$SLIDE)) 
	do
		python selection_seq_interval_bis.py $Results/Cortex_"$current".fasta $SCAFFOLD $i $j > "$current"_"$i"_"$j".fasta
		RES=$(minimap2 -x sr  $Results/Cortex_"$next".fasta "$current"_"$i"_"$j".fasta|head -n 1|cut -f 1-12)
		echo -e $i"\t"$j"\t"$RES |perl -pe 's/ /\t/g' >> mapping.txt
		rm "$current"_"$i"_"$j".fasta
	done

	# Combine the results
	mkdir -p $Results/minimap2/sliding_windows_mapping
	cat mapping.txt|awk 'NF > 2' > $Results/minimap2/sliding_windows_mapping/mapping_minimap2_sliding_windows_"$current"_"$next".txt
	rm mapping.txt

	#Plot the results
	Rscript ./Plotting_homology_new.R  $Results/minimap2/sliding_windows_mapping/mapping_minimap2_sliding_windows_"$current"_"$next".txt $current $next
	# Remove pointless files
	#rm minimap_plot.txt
	mv *pdf $Results/minimap2/sliding_windows_mapping/
	echo -E "END SLIDING WINDOWS $current $next"
done

echo -E "ALL SLIDING WINDOWS FINISHED"
