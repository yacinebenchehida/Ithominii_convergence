#!/bin/bash

##########################
# load necessary modules #
##########################
module load R/4.2.1-foss-2022a
module load Biopython/1.81-foss-2022b
module load MUMmer/3.23-GCCcore-9.3.0

##########################
# Set useful directories #
##########################
Results="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Results"
Inputs="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Inputs"
annotation="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Inputs/Optix_annotation.txt"

############################################################################################
# Find top SNPs in the GWAS and create a fasta in the peak regions and surrounding genes   #
############################################################################################
for i in $@
do
	ref_genome="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$i"
	fasta_sequence=$(ls $ref_genome|grep -E "*.fa$|*.fasta$")
	Scaffold=$(cat $annotation|grep $i|head -n 1|awk -F"\t" '{print $2}')
	starting_position=$(python3 ./find_min_max.py $annotation $i|awk '{print $1}')
	ending_position=$(python3 ./find_min_max.py $annotation $i|awk '{print $2}')
	echo $Scaffold $starting_position $ending_position
	python selection_seq_interval.py  $ref_genome/$fasta_sequence $Scaffold $starting_position $ending_position $i> $Results/Cortex_"$i".fasta
	#rm top_peak_positions_"$i".txt
done

echo FOUND PEAKS


######## Nucmer in sliding windows of 1kb
values=("$@")

for ((k = 0; k < ${#values[@]} - 1; k += 1)); do # Iterate through adjacent pairsS
	# Define the pair of species that are going to be analyses in this loop
	current="${values[k]}"
        next="${values[k + 1]}"
	echo -e "START SLIDING WINDOWS $current $next"
	echo -e "start\tend\tsubject\tsubjectLen\tsubjectStart\tsubjectEnd\tquery\tqueryLen\tqueryStart\tqueryEnd\tIdentity" >  mapping.txt
	# Define the size of the size genomic interval and the scaffold name that is going to be plotted
	SIZE=$(python3 -W ignore scaffold_size.py $Results/Cortex_"$current".fasta|awk '{print $2}')
	SCAFFOLD=$(python3 -W ignore scaffold_size.py $Results/Cortex_"$current".fasta|awk '{print $1}')
	# Define the window slide and the slide (here without overlap so both have the same size)
	WINDOWS=1000
	SLIDE=1000
	echo $SIZE $current

 	# Loop for each pair of species in the define scaffold over windows of the same size
	for ((i = 1, j = $WINDOWS; i < $SIZE && j < $SIZE; i = j + 1, j=j+$SLIDE)) 
	do
		echo -E "$i $j $SIZE"
		python selection_seq_interval_bis.py $Results/Cortex_"$current".fasta $SCAFFOLD $i $j > "$current"_"$i"_"$j".fasta
		nucmer --mum -c 20 -b 500 -l 10 --maxgap 500 -p tmp_"$i"_"$j" $Results/Cortex_"$next".fasta "$current"_"$i"_"$j".fasta
		show-coords -rcl tmp_"$i"_"$j".delta > tmp_nucmer_"$current"_"$next"_"$i"_"$j".txt
		(cat tmp_nucmer_"$current"_"$next"_"$i"_"$j".txt |grep -v "====="|awk 'NR> 4'|perl -pe 's/ +/\t/g' |perl -pe 's/^\t//g'|perl -pe 's/\|\t//g'|awk '{print $12"\t"$8"\t"$1"\t"$2"\t"$13"\t"$9"\t"$3"\t"$4"\t"$7}') > tmp3
		awk -v i="$i" -v j="$j" '{print i "\t" j "\t" $0}' tmp3 >> mapping.txt
		rm "$current"_"$i"_"$j".fasta tmp*
	done

	# Last windows
	j=$SIZE
	echo -E "$i $j $SIZE"
	python selection_seq_interval_bis.py $Results/Cortex_"$current".fasta $SCAFFOLD $i $j > "$current"_"$i"_"$j".fasta
        nucmer --mum -c 20 -b 500 -l 10 --maxgap 500 -p tmp_"$i"_"$j" $Results/Cortex_"$next".fasta "$current"_"$i"_"$j".fasta
       	show-coords -rcl tmp_"$i"_"$j".delta > tmp_nucmer_"$current"_"$next"_"$i"_"$j".txt
        (cat tmp_nucmer_"$current"_"$next"_"$i"_"$j".txt |grep -v "====="|awk 'NR> 4'|perl -pe 's/ +/\t/g' |perl -pe 's/^\t//g'|perl -pe 's/\|\t//g'|awk '{print $12"\t"$8"\t"$1"\t"$2"\t"$13"\t"$9"\t"$3"\t"$4"\t"$7}') > tmp3
        awk -v i="$i" -v j="$j" '{print i "\t" j "\t" $0}' tmp3 >> mapping.txt
        rm "$current"_"$i"_"$j".fasta tmp*

	# Combine the results
	mkdir -p $Results/mummer/sliding_windows_mapping
	cat mapping.txt|awk 'NF > 2' > $Results/mummer/sliding_windows_mapping/mapping_nucmer_sliding_windows_"$current"_"$next".txt
	rm mapping.txt

	#Plot the results
	Rscript ./plotting_nucmer_windows.R  $Results/mummer/sliding_windows_mapping/mapping_nucmer_sliding_windows_"$current"_"$next".txt $current $next $annotation
	# Remove pointless files
	#rm minimap_plot.txt
	mv *pdf $Results/mummer/sliding_windows_mapping/
	echo -E "END SLIDING WINDOWS NUCMER $current $next"
done

echo -E "ALL SLIDING WINDOWS FINISHED"
