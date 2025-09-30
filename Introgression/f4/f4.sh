#!/bin/bash
#SBATCH --time=00-01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=f4_f3

################
# Load modules #
################
module load BCFtools/1.19-GCC-13.2.0
module load  R/4.2.1-foss-2022a
module load PLINK/1.90-beta-7.4-x86_64

################
# Usage message
################
usage() {
    echo "Usage: $0 -v <vcf> -o <outdir> -p <pop_file> --p1 <pop1> --p2 <pop2> --p3 <pop3> [--p4 <pop4>] --stat <f3|f4>"
    echo ""
    echo "OPTIONS:"
    echo "  -v, --vcf      Input VCF file (bgzipped and indexed)."
    echo "  -o, --outdir   Output directory."
    echo "  -p, --pop      Population file (2 columns: popID sampleID)."
    echo "  --p1           Population name for P1."
    echo "  --p2           Population name for P2."
    echo "  --p3           Population name for P3."
    echo "  --p4           Population name for P4 (required only for f4)." Ignored if for f3. 
    echo "  --stat         Statistic to compute: 'f3' or 'f4'."
    echo "  -h             Show this help message and exit."
    exit 1
}

##################################
# Parse command-line arguments
##################################
if ! OPTIONS=$(getopt -o v:o:p:h --long vcf:,outdir:,pop:,p1:,p2:,p3:,p4:,stat: -n "$0" -- "$@"); then
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
        --stat) STAT="$2"; shift 2 ;;
        -h) usage ;;
        --) shift; break ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

##################################
# Check mandatory parameters
##################################
if [[ -z $VCF || -z $OUTDIR || -z $POP_FILE || -z $P1 || -z $P2 || -z $P3 || -z $STAT ]]; then
    echo "ERROR: Missing required arguments!"
    usage
fi

# Set P4 to NA if not provided (for f3 tests)
if [[ -z $P4 ]]; then
    P4="NA"
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
echo -e Stat:"${STAT}"
echo Options parsed

#################
# Define prefix #
#################
if [ "$STAT" == "f4" ]; then
    PREFIX=${P1}_${P2}_${P3}_${P4}
elif [ "$STAT" == "f3" ]; then
    PREFIX=${P1}_${P2}_${P3}
else
    echo "ERROR: STAT must be either 'f3' or 'f4'"
    exit 1
fi

RESULTS="$OUTDIR/$PREFIX"
echo -e The results are saved in: "${RESULTS}"
mkdir -p $RESULTS

#########################
# Create phenotype file #
#########################
if [ -f  $RESULTS/"$PREFIX"_phenotype_file.txt ]; then
    rm  $RESULTS/"$PREFIX"_phenotype_file.txt
fi
if [ -f  $RESULTS/"$PREFIX"_subspecies.txt ]; then
    rm  $RESULTS/"$PREFIX"_subspecies.txt
fi

touch $RESULTS/"$PREFIX"_phenotype_file.txt
touch $RESULTS/"$PREFIX"_subspecies.txt

for i in $P1 $P2 $P3 $P4; do
    if [ "$i" != "NA" ]; then
        awk -v pattern="$i" '$1 ~ pattern' $POP_FILE >> $RESULTS/"$PREFIX"_phenotype_file.txt
        echo $i >> $RESULTS/"$PREFIX"_subspecies.txt
    fi
done

####################################################
# Create subsetted VCF with the populations of interest
####################################################
bcftools view -S <(awk '{print $2}' "$RESULTS/${PREFIX}_phenotype_file.txt") \
  -i 'F_MISSING<0.2 && FORMAT/GQ>=10' \
  -m2 -M2 -v snps \
  -O z -o "$RESULTS/${PREFIX}.vcf.gz" \
  --threads 8 "$VCF"
tabix "$RESULTS/${PREFIX}.vcf.gz"
VCF="$RESULTS/${PREFIX}.vcf.gz"
echo Subsetted VCF generated

###############################################
# Run plink (needed for admixtools) #
###############################################
echo Running plink
plink --vcf "$VCF" --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --make-bed --out "$RESULTS/subspecies" --threads 8

echo Change contig names to numbers
awk '{if(!($1 in seen)) seen[$1]=++counter; $1=seen[$1]; print}' OFS='\t' $RESULTS/subspecies.bim > $RESULTS/Inputs.bim
mv $RESULTS/subspecies.bed $RESULTS/Inputs.bed
paste <(cat $RESULTS/"$PREFIX"_phenotype_file.txt|awk '{print $1}') <(awk '{$1=""; sub(/^ /,""); print}' $RESULTS/subspecies.fam) > $RESULTS/Inputs.fam
echo Plink done

##########################
# Run f3 or f4 statistics #
##########################
# Nota Bene: For f3, it will run f3(pop1; pop2, pop3), assuming one of the following interpretations:
# Option1: Pop1 is the "admixed" population and pop2 and pop3 are the source populations. In this case, a negative f3 indicates admixture.
# Option2: Pop1 is the outgroup population, and pop2 and pop3 are the populations being tested for genetic closeness.(then f3 is negative) or the outgroup ()
echo -e "Start $STAT using admixtools"
Rscript ./f4_f3.R $RESULTS/"$PREFIX"_subspecies.txt $RESULTS $STAT
echo DONE
