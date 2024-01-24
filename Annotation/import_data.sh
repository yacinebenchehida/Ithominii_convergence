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
SAMPLES=$($EDIRECT/esearch -db sra -query PRJNA725991 | $EDIRECT/efetch -format runinfo | cut -d "," -f 1,5,8,29|perl -pe 's/,/\t/g'|awk 'NR > 1 {print $1}')
echo $SAMPLES

for i in $SAMPLES; do echo $i; fastq-dump --split-files $i; done
for i in *fastq; do fastq-dump --split-3 $i; done
for i in RNA*; do pigz -d -p 12 $i; done
