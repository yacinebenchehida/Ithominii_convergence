#!/bin/bash
#SBATCH --time=00-01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=35G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=f4

################
# Load modules #
################
module load BCFtools/1.19-GCC-13.2.0
module load  R/4.2.1-foss-2022a
module load PLINK/1.90-beta-7.4-x86_64

################
# Script usage #
################
usage()
{
cat << EOF
Usage: ./f4.sh -v <vcf> -o <outdir> -p <pop_file> --p1 <pop1> --p2 <pop2> --p3 <pop3> --p4 <pop4> -h

OPTIONS:
  -v, --vcf        Input VCF file name (full dataset to be subsetted). The VCF should be bgzipped and tabixed.
  -o, --outdir     Output directory where results will be saved
  -p, --pop        Tab separate population file containing population (col1) and sample ID (exactly as in the VCF)
  --p1             Name of the first population (used as P1 in the f4 test)
  --p2             Name of the second population (used as P2 in the f4 test)
  --p3             Name of the third population (used as P3 in the f4 test)
  --p4             Name of the fourth population (used as P4 in the f4 test)
  -h               Show this usage information and exit
EOF
exit 0
}

##########################
# Parse input parameters #
##########################
OPTIONS=$(getopt -o v:o:p:h \
    --long vcf:,outdir:,pop:,p1:,p2:,p3:,p4: \
    -n "$0" -- "$@")

if [ $? -ne 0 ]; then
    usage
fi

eval set -- "$OPTIONS"

while true; do
    case "$1" in
        -v|--vcf) VCF="$2"; shift 2 ;;
        -o|--outdir) OUTDIR="$2"; shift 2 ;;
        -p|--pop) POP_FILE="$2"; shift 2 ;;
        --p1) P1="$2"; shift 2 ;;
        --p2) P2="$2"; shift 2 ;;
        --p3) P3="$2"; shift 2 ;;
        --p4) P4="$2"; shift 2 ;;
        -h) usage; shift ;;
        --) shift; break ;;
        *) break ;;
    esac
done

# Check mandatory parameters
if [[ -z $VCF || -z $OUTDIR || -z $POP_FILE || -z $P1 || -z $P2 || -z $P3 || -z $P4 ]]; then
    usage
fi

#################
# Define prefix #
#################
PREFIX=${P1}_${P2}_${P3}_${P4}
echo -e VCF file: "${VCF}"
echo -e Pop file: "${POP_FILE}"
echo -e P1:"${P1}"
echo -e P2:"${P2}"
echo -e P3:"${P3}"
echo -e P4:"${P4}"
echo Options parsed

#########################
# Create phenotype file #
#########################
RESULTS="$OUTDIR/$PREFIX"
mkdir -p $RESULTS
echo -e The results are saved in: "${RESULTS}"

if [ -f  $RESULTS/"$PREFIX"_phenotype_file.txt ]; then
    rm  $RESULTS/"$PREFIX"_phenotype_file.txt
fi
touch  $RESULTS/"$PREFIX"_phenotype_file.txt

for i in $P1 $P2 $P3 $P4; do cat ${POP_FILE} | awk -v pattern="$i" '$1 ~ pattern' >> $RESULTS/"$PREFIX"_phenotype_file.txt; done
for i in $P1 $P2 $P3 $P4; do echo $i >> $RESULTS/"$PREFIX"_subspecies.txt; done
echo POP FILE READY

####################################################
# Create subsetted vcf with the 4 taxa of interest #
####################################################
echo Subsetting VCF 
bcftools view -S <(awk '{print $2}' "$RESULTS/${PREFIX}_phenotype_file.txt") \
  -i 'F_MISSING<0.2 && FORMAT/GQ>=10' \
  -m2 -M2 -v snps \
  -O z -o "$RESULTS/${PREFIX}.vcf.gz" \
  --threads 8 "$VCF"
  
tabix "$RESULTS/${PREFIX}.vcf.gz"
echo Subsetted VCF generated

###############################################
# Run plink (output needed to run admixtools) #
###############################################
echo Running plink
plink --vcf "$VCF" --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --make-bed --out "$RESULTS/subspecies" --threads 8
echo Change contig names to numbers
awk '{if(!($1 in seen)) seen[$1]=++counter; $1=seen[$1]; print}' OFS='\t' $RESULTS/subspecies.bim > $RESULTS/Inputs.bim
mv $RESULTS/subspecies.bed $RESULTS/Inputs.bed
echo Plink done

##########################
# Compute f4 statistics #
#########################
echo Start f4 statistics using admixtools
Rscript ./f4.R $RESULTS/"$PREFIX"_subspecies.txt $RESULTS
echo DONE
