#!/bin/bash
#SBATCH --time=00-01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=35G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=f4

module load BCFtools/1.19-GCC-13.2.0

VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/multisp/Melinaea_menophilus/snps_NO_setGT_GQ10.vcf.gz"
POP_FILE="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/admixtools/f4/Inputs/species_samples.txt"
P1=$1
P2=$2
P3=$3
P4=$4
PREFIX=${P1}_${P2}_${P3}_${P4}

echo -e P1:"${P1}"
echo -e P2:"${P2}"
echo -e P3:"${P3}"
echo -e P4:"${P4}"

# Create phenotype file
RESULTS="../Results/$PREFIX"
mkdir -p $RESULTS

if [ -f  $RESULTS/"$PREFIX"_phenotype_file.txt ]; then
    rm  $RESULTS/"$PREFIX"_phenotype_file.txt
fi
touch  $RESULTS/"$PREFIX"_phenotype_file.txt

for i in $P1 $P2 $P3 $P4; do cat ${POP_FILE} |awk -v pattern="$i" '$1 ~ pattern'  >>  $RESULTS/"$PREFIX"_phenotype_file.txt; done
for i in $P1 $P2 $P3 $P4; do echo $i >> $RESULTS/"$PREFIX"_subspecies.txt; done
echo POP FILE READY

# Create subsetted vcf file for taxa that will be analysed
time bcftools view -S <(cat $RESULTS/"$PREFIX"_phenotype_file.txt|awk '{print $1}') $VCF| bcftools view -i 'F_MISSING < 0.2'| bcftools filter -e 'FORMAT/GQ < 10'| bcftools view -m2 -M2 -v snps -O z -o  $RESULTS/"$PREFIX".vcf.gz --threads 8
tabix $RESULTS/"$PREFIX".vcf.gz

# load plink
PLINK="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/plink_linux_x86_64_20231018/plink"

# Run plink to generate a *bim *fam *bed files used as input by admitools2
$PLINK --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --make-bed --out ../Inputs/subspecies 
awk '{if(!($1 in seen)) seen[$1]=++counter; $1=seen[$1]; print}' OFS='\t' ../Inputs/subspecies.bim > ../Inputs/Inputs.bim
mv ../Inputs/subspecies.bed ../Inputs/Inputs.bed

# run the f4 using admixtools
module load  R/4.2.1-foss-2022a
Rscript ./f4.R $RESULTS/"$PREFIX"_subspecies.txt $RESULTS
