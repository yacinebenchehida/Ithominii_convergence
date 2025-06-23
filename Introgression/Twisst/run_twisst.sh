#!/bin/sh
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=twisst

echo COMBINING ALL TREES INTO A SINGLE FILE
SCRIPTS=$9

###################
# Combine results #
###################
cd $1

if [ -f  topologies.txt ]; then
    rm   topologies.txt
fi
touch topologies.txt

if [ -f  missing_parts.txt ]; then
    rm   missing_parts.txt
fi
touch missing_parts.txt


cat "$2".txt | awk '{print $1}' | while read line; do
    dir="$2"_part_"$line"
    
    if [ -d "$dir" ]; then
        if ls "$dir"/* 1> /dev/null 2>&1; then
            cat "$dir"/* >> topologies.txt
        else
            echo "$line" >> missing_parts.txt
        fi
    else
        echo "$line" >> missing_parts.txt 
    fi
done

##################################################
# Remove lines genomics windows without topology #
##################################################
# Check if missing_parts.txt is empty
if [[ ! -s missing_parts.txt ]]; then
    # If missing_parts.txt is empty, output test_depenencies.txt unchanged
    awk '{print $2"\t"$3}' "$2".txt  > positions.txt
else
    # If missing_parts.txt is not empty, proceed with the filtering
    awk 'NR==FNR {remove[$1]; next} !($1 in remove)' missing_parts.txt "$2".txt  | awk '{print $2"\t"$3}' > positions.txt
fi

##############
# Run TWISST #
##############
# 1) Create G_OPTION
G_OPTION=$(cat *phenotype*| awk '{print $2}'|sort|uniq|perl -pe 's/(.+)\n/-g $1 /g')
echo $G_OPTION
echo G_OPTION CREATED

# 2) RUN TWISST
module load Biopython/1.83-foss-2023a
TWISST="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/twisst/twisst.py"
python $TWISST -t topologies.txt --method complete --groupsFile *phenotype* $G_OPTION --outputTopos Phylogeny.topos | gzip > Phylogeny.weights.tsv.gz
echo TWISST FINISHED

# 3) Plot TWISST RESULTS
module purge
module load  R/4.2.1-foss-2022a
Rscript $SCRIPTS/twisst.R $3 $4 0.012 $5 $6 $7 $8
echo TWISST RESULTS PLOTTED

############
# Cleaning #
############
rm -r "$2"_part_* *sh
echo TEMPORARY FILES CLEANED.
echo DONE
