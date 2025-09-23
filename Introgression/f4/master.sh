#!/bin/bash

# Define variables
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/multisp/Melinaea_menophilus/snps_NO_setGT_GQ10.vcf.gz" # Must be bgzipped and indexed  
OUTPUT="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/admixtools/f4/Results"
POP_FILE="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/admixtools/f4/Inputs/species_samples.txt" # Tab separate two columns text file (column 1: population; column 2: Sample_ID)

# # Arrays defining groups of four populations (P1, P2, P3, P4), where each index across the arrays specifies one f4 test comparison
P1_list=(pop11 pop12 pop13 pop14 pop15)
P2_list=(pop21 pop22 pop23 pop24 pop25)
P3_list=(pop31 pop32 pop33 pop34 pop35)
P4_list=(pop41 pop42 pop43 pop44 pop45)

# Loop over the indices of the population arrays; for each index, assign P1–P4 from the corresponding array entries and submit an f4 test job with those four populations
for i in ${!P1_list[@]}; do
    P1=${P1_list[$i]}
    P2=${P2_list[$i]}
    P3=${P3_list[$i]}
    P4=${P4_list[$i]}
    
    sbatch ./f4.sh -v $VCF -o $OUTPUT_PATH -p $POP_FILE --p1 ${P1} --p2 ${P2} --p3 ${P3} --p4 ${P4}   
done

