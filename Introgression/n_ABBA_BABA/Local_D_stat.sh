VCF="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/Twisst/6_Combine_intervals/Results/Melinaea_marsaeus/New/Melinaea_marsaeus.GQ20.DP5.snps.vcf.gz"
Input="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/D_statistics/Data"
Path_Output="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/D_statistics/Results"
Prefix="$1"
Chromsome="$2"
Start="$3"
End="$4"

module load bio/BCFtools/1.15.1-GCC-11.3.0


bcftools view --regions $Chromsome:$Start-$End $VCF | bcftools view -i 'F_MISSING < 0.30' > $Path_Output/"$Prefix".vcf
bgzip $Path_Output/"$Prefix".vcf
tabix $Path_Output/"$Prefix".vcf.gz


bcftools query -f '%POS  [ %GT]\n'] $Path_Output/"$Prefix".vcf.gz > $Path_Output/"$Prefix".genotypes.txt 

Rscript Local_D_stat.R /shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/D_statistics/ $Prefix menophilus_hicetas  menophilus_menophilus satevis_flavomacula lilis_imitata $Start $End
