#########################################################################################
# Split reference scaffolds into 99 pieces (or the number of scaffolds if less than 99) #
########################################################################################
./Split_ref_contigs.sh $1

#####################################
# Compute the size of each scaffold #
#####################################
module load bio/Biopython/1.79-foss-2021a
REF="/mnt/lustre/groups/biol-specgen-2018/kanchon/reference_genomes/$1"
FASTA_REF=$(ls /mnt/lustre/groups/biol-specgen-2018/kanchon/reference_genomes/$1|grep -E "*.fa$|*.fasta$")
DATA="/mnt/lustre/groups/biol-specgen-2018/yacine/Fixed_mutations_sliding_windows/Inputs/$1"

python -W ignore scaffold_size.py $REF/$FASTA_REF > $DATA/$1_contigs_size.txt

######################################################################
# Create a variable with the species name to create the the job names #
######################################################################
NUM_ARRAY=$(ls /mnt/lustre/groups/biol-specgen-2018/yacine/Split_reference_genomes/Results/$1/*part*|wc -l)
MAX_ARRAY=$((NUM_ARRAY - 1))
NAME=$(echo $1|cut -d "_" -f 2| cut -c 1-3)

###########################################################
# Compute the fixed mutation for each interval separately #
###########################################################
sbatch --array=0-$MAX_ARRAY --job-name=Fix_"$NAME" ./Fixed_sliding_windows.sh $@

###################################
# Extract the running job numbers #
###################################
running_jobs1=$(squeue|grep ybc502| grep Fix_"$NAME"| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')

#######################################################################           
# Combine the results in a single file and plot them along the genome #
#######################################################################
window_size=$(echo $@|perl -pe 's/([a-zA-Z_]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)(\s+)([0-9,]+)/$11/g')
sbatch --job-name=Combiner --dependency=aftercorr:$running_jobs1 ./Results_combiner.sh $1 $window_size



#Usage example with hydra vs histrio Chetone: ./master_script.sh Chetone_histrio 1,3,4,12,16,17,18,19,20,29,32,34,35,36,37 2,5,6,7,8,9,10,11,13,14,15,21,22,23,24,25,26,27,28,30,31,33,38,39,40,41,42 8 15 50000
