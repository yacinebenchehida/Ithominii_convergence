#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=twisst

RESULTS=$1
P1=$2
P2=$3
P3=$4
P4=$5
PO=$6
POP=$7
METHOD=$8

#######################################    
# 1 - Combine results in a single file #
########################################  
if [ "$METHOD" == "ML" ]; then
	cat $RESULTS/*/*.raxml.bestTree > $RESULTS/All_trees.txt
else
	cat $RESULTS/*/Rooted_NJ_tree.newick > $RESULTS/All_trees.txt
fi

############################    
# 2 - Prepare position file #
#############################
echo -e start"\t"end > $RESULTS/positions.txt


for ((i = 25956241, j = 25958240; i < 26031063 && j < 26031063; i = j + 1, j=j+2000))  ; do
	if [ "$METHOD" == "ML" ]; then 
		if ls $RESULTS/*"$i"*"$j"*/*bestTree 1> /dev/null 2>&1; then
        		echo -e $i"\t"$j >> $RESULTS/positions.txt
		fi
	elif [ "$METHOD" == "NJ" ]; then
		if ls $RESULTS/*"$i"*"$j"*/*newick 1> /dev/null 2>&1; then
			echo -e $i"\t"$j >> $RESULTS/positions.txt
		fi
	else
    		echo "Unknown method"
	fi
done

##################
# 3 - Run Twisst #
##################
TWISST="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/twisst/twisst.py"

touch  $RESULTS/groups.txt
for i in $P1 $P2 $P3 $P4 $PO; do cat ${POP}| grep $i >> $RESULTS/groups.txt; done
echo POP FILE READY


G_OPTION=$(cat $RESULTS/groups.txt| awk '{print $2}'|sort|uniq|perl -pe 's/(.+)\n/-g $1 /g')
python $TWISST -t $RESULTS/All_trees.txt --method complete --groupsFile $RESULTS/groups.txt $G_OPTION --outgroup $PO --outputTopos $RESULTS/Phylogeny.topos | gzip > $RESULTS/Phylogeny.weights.tsv.gz 
echo "Done running twisst"

###############################
# 4 - Plotting twisst results #
###############################
cp *.R  $RESULTS 
module load R/4.2.1-foss-2022a
Rscript --vanilla $RESULTS/twisst.R $RESULTS
echo "Done plotting twisst results"
