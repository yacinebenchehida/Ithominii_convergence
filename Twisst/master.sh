VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/multisp/Melinaea_marsaeus/Marsaeus_no_invariants_all_variations.vcf.gz"
GROUP="../Data/group_$1.txt"
CHR="SUPER_2"
P1="mothone"
P2="satevis"
P3="marsaeus_rileyi"
P4="marsaeus_phasiana"
PO="lilis"
METHOD="NJ"

############################################
# 1 - Run RAXML on each window of interest #
############################################
for ((i = 25956241, j = 25958240; i < 26031063 && j < 26031063; i = j + 1, j=j+2000))  ; do \
#for ((i = 23956241, j = 23958240; i < 26031063 && j < 26031063; i = j + 1, j=j+2000))  ; do \
sbatch ./Get_phylogenies.sh \
-v $VCF \
-c $CHR \
-s $i -e $j -f $GROUP \
-p1 $P1 -p2 $P2 -p3 $P3 -p4 $P4 -po $PO \
-m $METHOD \
-o ../Results/$1 -n "$CHR"_"$i"_"$j"_"$1"; \
done

#########################################
# 2 - summarize the results with twisst #
#########################################
running_jobs1=$(squeue|grep ybc502| grep raxml| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
sbatch  --dependency=aftercorr:$running_jobs1  ./twisst.sh  ../Results/$1 $P1 $P2 $P3 $P4 $PO $GROUP $METHOD
