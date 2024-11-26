VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/multisp/Melinaea_menophilus/Multisp_menophilus.DP4_GQ5_QUAL5_invariants.vcf.gz"

# Cortex
sbatch ./MuteBass.sh -v $VCF \
-c SUPER_18 -s 1 -e 16520369 \
-p /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/MuteBass/Inputs/yellow_band_species.txt \
-t "(1,((2,(3,4)),5))" \
-f 0.5 \
-o /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/MuteBass/Results/SUPER_18

# Optix

# Antennapedia
