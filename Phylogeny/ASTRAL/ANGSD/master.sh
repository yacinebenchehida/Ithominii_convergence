#!/bin/bash

#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=ast_mas
#SBATCH --time=0-02:00:00

##########################
# load necessary modules #
##########################
module load R/4.2.1-foss-2022a
module load Biopython/1.81-foss-2022b
module load HTSlib/1.15.1-GCC-11.3.0
module load Java/13.0.2

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
# Split fasta files into windows of 100000 kb #
###############################################

######################################
# Get ML or NJ tree for each windows #
######################################

############################################
# Combine all the trees into a single file #
############################################

##############
# Run ASTRAL #
##############
# Forcing monophyly of known species

# Without forcing the monophyly


####################
# Plot ASTRAL tree #
####################
