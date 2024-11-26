#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=MutBa

#################
# Load packages #
#################
module load BCFtools/1.19-GCC-13.2.0
module load VCFtools/0.1.16-GCC-11.2.0
module load Pysam/0.22.0-GCC-12.3.0

###################################
# Set usage and script parameters #
###################################
usage()
{
cat << EOF
./test_mutebass.sh -v <vcf> -c <chromosome> -s <star_position> -e <end_position> -p <pop_file> -t <tree> -f <target_frequency> -o <path_output> -h

OPTIONS:
  -v            Input VCF file name
  -c		Contig/Scaffold/Chromosome name
  -s		Starting position in the VCF
  -e		End position in the VCF
  -p            Population file
  -t		List of taxa to analyses (comma-separated)
  -f		Target frequency for testing balancing selection
  -o		Output folder  path (absolute file only)
  -h            usage information and help (this message)
EOF
exit 0
}

# Parse command-line options
options=$(getopt -o v:c:s:e:p:t:f:o:h -- "$@")

# Check for getopt errors
if [ $? -ne 0 ]; then
  usage
fi

# Set the parsed options back to positional parameters
eval set -- "$options"

# Extract options and their arguments
while [[ $# -gt 0 ]]
do
    case $1 in
		-v)
		  VCF=$2
		  shift 2
		  ;;
		-c)
		  CHR=$2
		  shift 2
		  ;;
		-s)
		  START=$2
		  shift 2
		  ;;
		-e)
		  END=$2
		  shift 2
		  ;;
		-p)
		  POP=$2
		  shift 2
		  ;;
		-t)
		  TREE=$2
		  shift 2
		  ;;
                -f)
                  FREQ=$2
                  shift 2
                  ;;
		-o)
		  OUTPUT=$2
		  shift 2
		  ;;
		-h)
		  usage
		  ;;
		--)
		  shift
		  break
		  ;;
		*)
		  usage
		  ;;
	esac
done

# Check if all mandatory arguments are provided
if [[ -z $VCF ]] || [[ -z $CHR ]] || [[ -z $START ]] || [[ -z $END ]] || [[ -z $POP ]] || [[ -z $TREE ]] || [[ -z $FREQ ]] || [[ -z $OUTPUT ]] 
then
	usage
	exit 1
fi

#########################
# 2 - Display arguments #
#########################
echo -e VCF file:"${VCF}"
echo -e Chromosome:"${CHR}"
echo -e Start:"${START}"
echo -e End:"${END}"
echo -e Population file:"${POP}"
echo -e TREE:"${TREE}"
echo -e Target Frequency: "${FREQ}"
echo -e OUTPUT:"${OUTPUT}"

############################
# Create working directory #
############################
if [ -d ${OUTPUT} ]; then
    rm -r ${OUTPUT}
fi
mkdir -p ${OUTPUT}

# Ingroup
bcftools view --regions $CHR:$START-$END -S <(cat $POP|awk '{print $1}') $VCF | bcftools view -U -i 'TYPE=="snp"' -m2 -M2 -v snps -O v|bcftools view -i 'F_MISSING < 0.5'  > $OUTPUT/test.vcf
bgzip $OUTPUT/test.vcf
tabix $OUTPUT/test.vcf.gz
echo INGROUP VCF READY

# outgroup
bcftools view --regions $CHR:$START-$END -T <(bcftools query -f '%CHROM\t%POS\n' $OUTPUT/test.vcf.gz) -S <(cat $POP|awk '$2=="outgroup" {print $1}') $VCF > $OUTPUT/outgroup.vcf
bgzip $OUTPUT/outgroup.vcf
tabix $OUTPUT/outgroup.vcf.gz
echo OUTGROUP VCF READY

# Allele frequencies outgroup
indv_command=$(for i in $(bcftools query -l  $OUTPUT/outgroup.vcf.gz); do echo -e "--indv $i"; done |perl -pe 's/\n/ /g')
vcftools --gzvcf $OUTPUT/outgroup.vcf.gz --freq $indv_command --out $OUTPUT/outgroups
echo ALLELE FREQUENCIES FOR OUTGROUP READY

# Set outgroup state as ancestral state (when data are available i.e. no missing data)
awk '{
    split($5, ref, ":"); 
    split($6, alt, ":"); 
    if (ref[2] == "-nan") { 
        print $2, "nan"; 
    } else if (ref[2] > 0.5) { 
        print $2, ref[1]; 
    } else if (ref[2] == 0.5) { 
        print $2, "nan"; 
    } else { 
	print $2, alt[1]; 
    }
}' <(grep -v CHROM $OUTPUT/outgroups.frq) > $OUTPUT/ancestral.alleles
echo ANCESTRAL ALLELES FROM OUTGROUP SET

# Select 5 random samples for each taxa
python3 select_random.py $POP $OUTPUT
echo RANDOM SAMPLES SELECTED

# Get allele frequencies for all subsamples (for each species)
for i in $(cat $POP|awk '{print $2}'|grep -v outgroup| sort -u)
do
	echo $i
	bcftools view --regions $CHR:$START-$END  -T <(bcftools query -f '%CHROM\t%POS\n' $OUTPUT/test.vcf.gz) -S <(cat $OUTPUT/"$i".txt|awk '{print $1}') $VCF > $OUTPUT/"$i".vcf
	bgzip $OUTPUT/"$i".vcf
	tabix $OUTPUT/"$i".vcf.gz
	vcftools --gzvcf $OUTPUT/"$i".vcf.gz --freq --out $OUTPUT/"$i".frequencies
	grep -f <(cat $OUTPUT/ancestral.alleles|grep nan |awk '{ print $1}') $OUTPUT/"$i".frequencies.frq >  $OUTPUT/"$i".frq
	bcftools query -f '%POS\t%REF\t%ALT\t[%GT\t]\n' $OUTPUT/"$i".vcf.gz > $OUTPUT/"$i".geno
	rm $OUTPUT/"$i".frequencies.frq
done
echo ALLELE FREQUENCIES FOR FOR SUBSET DEFINED

# Create MuteBass input
SPECIES=$(cat $POP|awk '{print $2}'|grep -v outgroup|sort -u|perl -pe 's/\n/,/g'|perl -pe 's/,$//g')
echo $SPECIES

module load R/4.2.1-foss-2022a
Rscript ./Input_MuteBass.r $SPECIES $OUTPUT
echo MUTEBASS RAW INPUT READY

# Remove Snps not compatible with the with the species phylogeny
python3 /mnt/scratch/projects/biol-specgen-2018/yacine/Tools/MuteBaSS/MuteBaSS.py -c 1,2 -i $OUTPUT/input_mutebass --tree $TREE --check
mv snp_2_remove.txt $OUTPUT
python3 remove_line_not_following_tree.py $OUTPUT/input_mutebass $OUTPUT/snp_2_remove.txt $OUTPUT/Final_mutebass_input.txt 
echo MUTEBASS CLEAND INPUT READY

# Get HKA configuration file
python3 /mnt/scratch/projects/biol-specgen-2018/yacine/Tools/MuteBaSS/MuteBaSS.py -i $OUTPUT/Final_mutebass_input.txt -c 1,2 --getConfig --config $OUTPUT/HKA_input 

#  RUN NCD and HKA
python3 /mnt/scratch/projects/biol-specgen-2018/yacine/Tools/MuteBaSS/MuteBaSS.py -c 1,2 -i $OUTPUT/Final_mutebass_input.txt --fixSize -w 10000 -s 1000 -o $OUTPUT/final_results.txt --HKA --config $OUTPUT/HKA_input --NCDopt --NCD --NCDsub --tf $FREQ
echo NCD DONE

#  RUN NCDmid
python3 /mnt/scratch/projects/biol-specgen-2018/yacine/Tools/MuteBaSS/MuteBaSS.py -c 1,2 -i $OUTPUT/Final_mutebass_input.txt --fixSize -w 10000 -s 1000 -o $OUTPUT/final_results.txt --NCDmid --tf $FREQ
echo NCD MID DONE

#  Plot results 
Rscript ./Plotting_results.r $OUTPUT/final_results.txt $OUTPUT/final_results_fixSizeCT.txt 
mv *pdf $OUTPUT
# USAGE EXAMPLE: ./test_mutebass.sh -v /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/multisp/Melinaea_menophilus/multisp_Melinaea_menophilus_genotypeGVCF.intervals_55.filters.DP4_GQ5_QUAL5.invariants.vcf.gz -c SUPER_18 -s 13000000 -e 13800000 -p $POP -t "(1,((2,(3,4)),5))" -o /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/MuteBass/Results/testons_la_vie_en_rose
