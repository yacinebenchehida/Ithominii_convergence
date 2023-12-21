#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=10G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=Comb_int

###########################
# 0 - DEFINE USEFUL PATHS #
###########################
INT_VCFS="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/$1"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/$1"

#############################   
# 1 - CREATE RESULTS FOLDER #
############################# 
mkdir -p $RESULTS 

#############################   
# 2 - COMBINE VCF INTERVALS #
############################# 
cat $INT_VCFS/*intervals_0.* > $RESULTS/"$1".GQ5.DP4_FM50.snps.vcf
for j in {1..59}; do echo interval "$j"; cat $INT_VCFS/*intervals_"$j".*|grep -v '#' >> $RESULTS/"$1".GQ5.DP4_FM50.snps.vcf; done

#############################
# 3 - COMPRESSED FINAL VCFS #
############################# 
module load  tabixpp/1.1.0-GCC-11.2.0
bgzip -@ 12 $RESULTS/"$1".GQ20.DP5.snps.vcf && tabix -p vcf $RESULTS/"$1".GQ20.DP5.snps.vcf.gz
