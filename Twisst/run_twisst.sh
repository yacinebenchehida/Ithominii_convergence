echo COMBINING ALL TREES INTO A SINGLE FILE

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
grep -v -F -f missing_parts.txt topologies.txt > my_final_topologies.txt
rm topologies.txt

##############
# Run TWISST #
##############
# 1) Create G_OPTION
G_OPTION=$(cat *phenotype*| awk '{print $2}'|sort|uniq|perl -pe 's/(.+)\n/-g $1 /g')
echo $G_OPTION
echo G_OPTION CREATED

# 2) RUN TWISST
TWISST="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/twisst/twisst.py"
python $TWISST -t my_final_topologies.txt --method complete --groupsFile *phenotype* $G_OPTION --outputTopos Phylogeny.topos | gzip > Phylogeny.weights.tsv.gz

# 3) Plot TWISST RESULTS
Rscript $SCRIPTS/twisst_zscore.R $3 $4 0.020 $5 $6 $7 $8
