#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G 
#SBATCH --account=BIOL-SPECGEN-2018 

################
# Define paths #
################
SCRIPT="/mnt/lustre/groups/biol-specgen-2018/kanchon/script_header_files/fixed_diffs_hom_get_from_vcf2.pl"
PATH_RESULTS="/mnt/lustre/groups/biol-specgen-2018/yacine/Fixed_mutations_sliding_windows/Results/$1"
REF="/mnt/lustre/groups/biol-specgen-2018/kanchon/reference_genomes/$1"
FASTA_REF=$(ls /mnt/lustre/groups/biol-specgen-2018/kanchon/reference_genomes/$1|grep -E "*.fa$|*.fasta$")
PATH_DATA="/mnt/lustre/groups/biol-specgen-2018/yacine/Fixed_mutations_sliding_windows/Inputs/$1"
#VCF="/users/ybc502/scratch/Conv_Evol/Bio_inf/6_Combine_intervals/Results/$1/*vcf.gz"
VCF="/users/ybc502/scratch/Conv_Evol/Analyses_bio/Twisst/6_Combine_intervals/Results/Melinaea_marsaeus/Melinaea_marsaeus.GQ20.DP5.snps.vcf.gz"
#VCF="/users/ybc502/scratch/Conv_Evol/Bio_inf/6_Combine_intervals/Results/All_marsaeus/Melinaea_marsaeus.GQ20.DP5.snps.vcf.gz"
CONT="/mnt/lustre/groups/biol-specgen-2018/yacine/Split_reference_genomes/Results/$1/contigs_$1_part_*0${SLURM_ARRAY_TASK_ID}"
CONTIGS=$(echo $CONT|awk '{print $1}')

##############################
# Create working directories #
##############################
mkdir -p $PATH_RESULTS
mkdir -p $PATH_DATA

####################
# Import libraries #
####################
module load bio/BCFtools/1.15-GCC-11.2.0 
module load bio/SAMtools/1.15-GCC-11.2.0
module load  bio/tabix/0.2.6-GCCcore-7.3.0

#######################################
# Define arguments for Kanchon script #
#######################################
arg1=$(echo $@|perl -pe 's/([a-zA-Z_]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)/$3/g')
arg2=$(echo $@|perl -pe 's/([a-zA-Z_]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)/$5/g')
arg3=$(echo $@|perl -pe 's/([a-zA-Z_]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)/$7/g')
arg4=$(echo $@|perl -pe 's/([a-zA-Z_]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)/$9/g')

window_size=$(echo $@|perl -pe 's/([a-zA-Z_]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)/$11/g')


echo "Group1: "$arg1 
echo "Group2: "$arg2 
echo "Minimum individuals Group1: "$arg3 
echo "Minimum individuals Group2: "$arg4
echo "window size: "$window_size

#################################################################################
# Create a tmp vcf files with a subset of all contigs (to speed up the analyses #
################################################################################# 
cat $CONTIGS|while read line # While read each contig in the interval
	do bcftools view --regions $line $VCF > $PATH_RESULTS/$line.vcf # Create a vcf for the contig of interest
	bgzip $PATH_RESULTS/$line.vcf
	tabix $PATH_RESULTS/$line.vcf.gz
	CONTIG_SIZE=$(cat /mnt/lustre/groups/biol-specgen-2018/yacine/Fixed_mutations_sliding_windows/Inputs/$1/$1_contigs_size.txt|grep -P "$line\t" |awk '{print $2}') # Get contig size
	for ((i = 1, j = $window_size; i < $CONTIG_SIZE && j < $CONTIG_SIZE; i = j + 1, j=j+$window_size)) # subset a window of n kb
		do bcftools view --regions "$line:$i-$j" $PATH_RESULTS/$line.vcf.gz > $PATH_RESULTS/"$line"_"$i"_"$j".vcf  # for each window create a vcf file for the corresponding window
		echo $i $j $line
		$SCRIPT $arg1 $arg2 $arg3 $arg4 $PATH_RESULTS/"$line"_"$i"_"$j".vcf # Calculate fixed mutation for each snp in the window
		fixed=$(cat $PATH_RESULTS/"$line"_"$i"_"$j".vcf_minind*diff|grep -v "#CHROM"|wc -l) # Sum up all the fixed mutations in the window
		het_hem=$(cat $PATH_RESULTS/"$line"_"$i"_"$j".vcf_minind*het|grep -v "#CHROM"|wc -l) # Sum up all the hetero vs homo  mutations in the window
		echo -e $line"\t"$i"\t"$j"\t"$fixed"\t"$het_hem >> $PATH_RESULTS/"$line"_fixed_mutations.txt # print the result to a file for that window
		rm $PATH_RESULTS/"$line"_"$i"_"$j"*
	done 
	j=$CONTIG_SIZE # Do the same for the last window
	echo $i $j $line
	bcftools view --regions "$line:$i-$j" $PATH_RESULTS/$line.vcf.gz > $PATH_RESULTS/"$line"_"$i"_"$j".vcf 
	$SCRIPT $arg1 $arg2 $arg3 $arg4 $PATH_RESULTS/"$line"_"$i"_"$j".vcf
	fixed=$(cat $PATH_RESULTS/"$line"_"$i"_"$j".vcf_minind*diff|grep -v "#CHROM"|wc -l)
	het_hem=$(cat $PATH_RESULTS/"$line"_"$i"_"$j".vcf_minind*het|grep -v "#CHROM"|wc -l) # Sum up all the hetero vs homo  mutations in the window 
	echo -e $line"\t"$i"\t"$j"\t"$fixed"\t"$het_hem >> $PATH_RESULTS/"$line"_fixed_mutations.txt
	rm $PATH_RESULTS/"$line"_"$i"_"$j"*
	rm $PATH_RESULTS/"$line".vcf*
 done
