#!/bin/bash
#SBATCH --job-name=busco        # Job name
#SBATCH --partition=nodes
#SBATCH --ntasks=1                       # Run on a single CPU
#SBATCH --cpus-per-task=8                # Number of CPU cores per task
#SBATCH --mem=30gb                        # Job memory request
#SBATCH --time=04:00:00                  # Time limit hrs:min:sec
#SBATCH --account=BIOL-SPECGEN-2018       # Project account

############################################
# 0 - DEFINE WORKING DIRECTORIES AND PATHS #
############################################
RESULTS="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/BUSCO/Results/$1"
ref_genome="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/BUSCO/Inputs/"
fasta_sequence=$(ls $ref_genome|grep -E "*.fa$|*.fasta$|*.fna$"|grep -E $1)
mkdir -p $RESULTS

##################    
# 1 - LOAD BUSCO #
##################    
module load bio/BUSCO/5.4.3-foss-2020b

#################       
# 2 - RUN BUSCO #
#################   
/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/busco/bin/busco -m genome -i $ref_genome/$fasta_sequence -o $1 --out_path $RESULTS -l lepidoptera_odb10 -c 1 --config  /shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/busco/config/config.ini
