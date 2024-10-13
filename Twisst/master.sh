VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/multisp/Melinaea_marsaeus/multisp_Melinaea_marsaeus_genotypeGVCF.intervals_8.filters.DP4_GQ5_QUAL5.invariants.vcf.gz"
#VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/multisp/Melinaea_marsaeus/Marsaeus_no_invariants_all_variations.vcf.gz"
GROUP="../Data/group_$1.txt"
CHR="SUPER_2"
START=24000000
END=27000000
TAXA="mothone,satevis,marsaeus_rileyi,marsaeus_phasiana"
OUTPUT_PATH="."
PREFIX="test"
POP="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/Relate/Inputs/hindwing_black.txt"
METHOD="ML"

###################################
# 2' - SOFTWARE LAODING AND PATHS #
###################################
module load BCFtools/1.19-GCC-13.2.0
module load R/4.2.1-foss-2022a
module load VCFtools/0.1.16-GCC-11.2.0
module load Biopython/1.83-foss-2023a
SHAPEIT="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/shapeit4/bin/shapeit4.2"
TWISST="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/twisst/twisst.py"
SCRIPTS=$(pwd)

##############################
# 3 - Get the phenotype file #
##############################
# Create working directory
if [ -d $OUTPUT_PATH/$PREFIX ]; then
    rm -r $OUTPUT_PATH/$PREFIX
fi
mkdir -p $OUTPUT_PATH/$PREFIX

# Create phenotype file
if [ -f  $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt ]; then
    rm  $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt
fi
touch  $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt

SPECIES=$(echo $TAXA|perl -pe 's/,/ /g')
echo $SPECIES
for i in $SPECIES; do cat ${POP}| grep -P "\s$i\s*(\w)*\s*$" >>  $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt; done
echo POP FILE READY


#######################################################################
# Create VCF for genomics regions of interest and species of interest #
#######################################################################
bcftools view --regions $CHR:$START-$END -S <(cat $OUTPUT_PATH/$PREFIX/"$PREFIX"_phenotype_file.txt|awk '{print $1}') $VCF > $OUTPUT_PATH/$PREFIX/"$PREFIX".vcf
bgzip $OUTPUT_PATH/$PREFIX/"$PREFIX".vcf
tabix $OUTPUT_PATH/$PREFIX/"$PREFIX".vcf.gz
echo "SUBSETTED VCF FILE READY"

# Create VCFs with 100 SNPs each
echo "STARTED SPLITTING VCF INTO SMALLER VCF WITH 100 SNPS EACH"
python3 ./splice.py $OUTPUT_PATH/$PREFIX/"$PREFIX".vcf.gz  $OUTPUT_PATH/$PREFIX/$PREFIX 100 $OUTPUT_PATH/$PREFIX/"$PREFIX".txt
echo "SPLITTING FINISHED"

# gzip and infer phylogeny for each vcf
python3 ./Phylogenies_job_submitter.py $PREFIX $METHOD $OUTPUT_PATH/$PREFIX $OUTPUT_PATH/$PREFIX
