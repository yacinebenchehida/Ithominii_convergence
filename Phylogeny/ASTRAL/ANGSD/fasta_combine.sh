#!/bin/bash

#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-00:15:00

module load Biopython/1.81-foss-2022b

RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/ASTRAL/Results"
cd $RESULTS
for i in *
do
	cd $i
	touch "$i"_combined.fasta
	for fasta in *gz
	do
		id=$(echo $fasta|perl -pe 's/_'"$i"'.fa.gz//g')
		python3 /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/ASTRAL/Scripts/Add_sample_id_to_header.py $fasta $id > "$id"_"$i".fas
		cat "$id"_"$i".fas >> "$i"_combined.fasta
		rm $fasta "$id"_"$i".fas
	done
	cd $RESULTS
done 
