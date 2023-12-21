#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G  
#SBATCH --job-name=fd


###############################
# 0 - Load necessry libraries #
###############################
module load bio/HTSlib/1.17-GCC-12.2.0
module load bio/BCFtools/1.15.1-GCC-11.3.0

#########################################
# 1 - DEFINE USEFUL PATHS AND VARIABLES #
#########################################
Prefix="$1"
POP_FILE="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/fd_statistics/Inputs/$2"
Path_Output="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/fd_statistics/Results/$Prefix"
SCRIPT="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/genomics_general"
VCF="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/Twisst/6_Combine_intervals/Results/Melinaea_marsaeus/New/Melinaea_marsaeus.GQ20.DP5.snps.vcf.gz"
Chromsome="$3"
Start="$4"
End="$5"
P1="$6"
P2="$7"
P3="$8"
PO="$9"
window="${10}"

mkdir -p $Path_Output

################################################
# 2 - GENERATE A VCF FOR THE REGION TO ANALYSE #
################################################   
bcftools view --regions $Chromsome:$Start-$End $VCF | bcftools view -i 'F_MISSING < 0.80' > $Path_Output/"$Prefix".vcf
bgzip $Path_Output/"$Prefix".vcf
tabix $Path_Output/"$Prefix".vcf.gz

#########################################################################                      
# 3 - Create the geno file necessary to get the sliding window analyses #
#########################################################################
python3 $SCRIPT/VCF_processing/parseVCF.py \
-i $Path_Output/"$Prefix".vcf.gz -o $Path_Output/"$Prefix".geno.gz
#python3 $SCRIPT/VCF_processing/parseVCF.py -i $VCF -o $Path_Output/"$Prefix".geno.gz

##############################       
# 4 - Get the phenotype file #
############################## 
for i in $P1 $P2 $P3 $PO; do cat $POP_FILE| grep $i >> $Path_Output/"$Prefix"_phenotype_file.txt; done

#############################
# 5 - fd in sliding windows #
#############################
python3 $SCRIPT/ABBABABAwindows.py \
-g $Path_Output/"$Prefix".geno.gz \
-f phased \
-o $Path_Output/output.csv \
--windType sites \
-w $window -m 1 \
-P1 $P1 -P2 $P2 -P3 $P3 -O $PO \
-T 4 \
--minData 0.2 \
--popsFile $Path_Output/"$Prefix"_phenotype_file.txt \
--writeFailedWindows

cat $Path_Output/output.csv | grep -v nan|perl -pe 's/,/\t/g' > $Path_Output/Results_"$Chromsome"_"$Start"_"$End"_"$P1"_"$P2"_"$P3"_"$PO".txt

#############################
# 6 - Plotting ABBA and fd  #
#############################
Rscript ./Statistics.R /shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/fd_statistics/  $Prefix Results_"$Chromsome"_"$Start"_"$End"_"$P1"_"$P2"_"$P3"_"$PO".txt $window $P1 $P2 $P3 $PO
 

# USAGE: ./fd.sh test3 Phenotype_Apex_spot.txt SUPER_5 20400000 20800000 marsaeus_phasiana marsaeus_rileyi satevis_maeolus lilis_parallelis 50
