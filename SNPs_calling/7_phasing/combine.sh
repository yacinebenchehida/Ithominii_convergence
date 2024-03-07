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
INT_VCFS="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/9_Phasing/Results/$1"

#############################   
# 2 - COMBINE VCF INTERVALS #
############################# 
zcat $INT_VCFS/*intervals_1.*snps.phased.vcf.gz > $INT_VCFS/"$1".DP4_GQ5_QUAL5_F20.snps.phased.vcf
for j in $(seq 2 $1); do echo interval "$j"; zcat $INT_VCFS/*intervals_"$j".*snps.phased.vcf.gz|grep -v '#' >> $INT_VCFS/"$1".DP4_GQ5_QUAL5_F20.snps.phased.vcf; done

#############################
# 3 - COMPRESSED FINAL VCFS #
############################# 
module load  tabixpp/1.1.0-GCC-11.2.0
bgzip -@ 12 $INT_VCFS/"$1".DP4_GQ5_QUAL5_F20.snps.phased.vcf && tabix -p vcf $INT_VCFS/"$1".DP4_GQ5_QUAL5_F20.snps.phased.vcf.gz
