#!/bin/bash

#SBATCH --mem=3GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=ast_mas
#SBATCH --time=0-00:50:00

##########################
# load necessary modules #
##########################
module load R/4.2.1-foss-2022a
module load Biopython/1.81-foss-2022b
module load HTSlib/1.15.1-GCC-11.3.0
module load Java/13.0.2
module load RAxML/8.2.12-gompi-2021a-hybrid-avx2
module load Armadillo/12.6.2-foss-2023a

#########################
# Get List of scaffolds #
#########################
python ./scaffold_size.py /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/Hypothyris_anastasia/hypothyris_anastasia_mtDNA_10_10_23.fasta |awk '{print $1}' > ../Inputs/list_scaffold.txt

###############################################
# For each scaffold convert fasta to bam file #
###############################################
sbatch ./bam2fasta.sh 

############################################################
# Combine fasta files of each scaffold in a different file #
############################################################
# Extract the running job numbers 
running_jobs1=$(squeue|grep ybc502| grep bam2fas| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
echo $(eval echo "$running_jobs1")
sbatch --job-name=CombFasta --dependency=aftercorr:$running_jobs1 ./fasta_combine.sh 

###############################################
# Split fasta files into windows of 100000 kb and get the ML tree for each window and root the trees
###############################################
# Extract the running job numbers 
running_jobs2=$(squeue|grep ybc502| grep CombFast| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
echo $(eval echo "$running_jobs2")
sbatch --job-name=Sliding --dependency=aftercorr:$running_jobs2 ./run_sliding_windows_raxml.sh 

######################################
# Run ASTRAL and plot the final tree #
######################################
# Forcing monophyly of known species

# Without forcing the monophyly
running_jobs3=$(squeue|grep ybc502| grep Sliding| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
echo $(eval echo "$running_jobs3")
sbatch --dependency=aftercorr:$running_jobs3 ./ASTRAL.sh
