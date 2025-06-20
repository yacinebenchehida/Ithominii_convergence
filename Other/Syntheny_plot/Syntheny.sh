################
# Load modules #
################
module purge
module load bio/Biopython/1.81-foss-2022b
module load bio/BLAST+/2.14.0-gompi-2022b
module load R/4.2.1-foss-2022a

######################################################################
# Set paths (THIS PART NEEDS TO BE MODIFIED IF USED IN ANOTHER USER) #
######################################################################
RESULTS="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/syntheny/Results"
INPUTS="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/syntheny/Inputs"
DB="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/NCBI/uniprot"

##################################################
# Go to the directory containing the fasta files #
##################################################
cd $INPUTS

######################################################
# Check if the results file exist and if so erase it #
######################################################
[ -e $RESULTS/results.txt ] && rm $RESULTS/results.txt

##################################
# Run the pipeline for each gene #
##################################
for gene in *fasta # Loop over gene fasta file
	do for sp in Chetone_histrio Heliconius_melpomene # Loop over species: PLEASE CHANGE MANUALLY THE SPECIES YOU WANT TO ANALYSE
		do 
		echo -e $gene"\t"$sp
		tblastn -query $gene -db $DB/"$sp" -num_threads 8 -outfmt 7 |grep -v "#"| sort -k 12 -n -r|head -n 1 |awk '{print $2"\t"$12"\t"$9"\t"$10}' # Blast each 
		MID_POINT=$(tblastn -query $gene -db $DB/"$sp" -num_threads 8 -outfmt 7 |grep -v "#"| sort -k 12 -n -r|head -n 1 |awk '{print ($10+$9)/2}') # Get mid points for the best local hit
		GENE=$(echo $gene|cut -d "." -f 1) # Get gene name
		python3  -W ignore ../Script/scaffold_size.py $INPUTS/$gene $MID_POINT > $RESULTS/gene_"$gene"_interval_"$sp".txt $GENE # Get roughly the position of the genes based on the position of the mid point
		paste <(echo -e $sp"\t"$GENE) <(cat $RESULTS/gene_"$gene"_interval_"$sp".txt) >> $RESULTS/results.txt # Combine the results in a file called results.txt
		rm $RESULTS/gene_"$gene"_interval_"$sp".txt # Remove intermediate files 
	done
done

############################
# Make syntheny plots in R #
############################
FIRST_GENE=$(cat $RESULTS/results.txt|head -n 1|awk '{print $2}')
Rscript ../Script/syntheny.R $RESULTS $FIRST_GENE
