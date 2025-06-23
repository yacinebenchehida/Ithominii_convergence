# Author: Yacine Ben Chehida

#SBATCH --time=08:01:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=10G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=imp_data

EDIRECT="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/edirect"
module load SRA-Toolkit/3.0.3-gompi-2022a
module load pigz/2.7-GCCcore-12.2.0

# Figure out samples to download
SAMPLES=$($EDIRECT/esearch -db sra -query PRJEB60199 | $EDIRECT/efetch -format runinfo | cut -d "," -f 1,5,8,29|perl -pe 's/,/\t/g'|awk 'NR > 1 {print $1}'|grep -E "ERR11837472|ERR12245539")
echo $SAMPLES

# Download sample and transform them into R1 and R2 reads
for i in $SAMPLES; do echo $i; fastq-dump --split-files $i; done

# Combine the R1 reads of all samples into a single file (same for R2).
for i in *_1.fastq; do echo $i; cat $i >> RNA_seq_Melinaea_R1.fastq ; done
for i in *_2.fastq; do echo $i; cat $i >> RNA_seq_Melinaea_R2.fastq ; done

# gzip the R1 and R2 fastq files
for i in RNA*; do pigz -p 12 $i; done
