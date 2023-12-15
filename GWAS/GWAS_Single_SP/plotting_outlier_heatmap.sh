SPECIES="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Inputs/$1/Grouping_"$2".txt"
PHENOTYPES="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Inputs/$1/GEMMA_encoding_phenotype_"$2".txt"
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/$1/*.vcf.gz"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Results/$1/$2"
Threshold=$3

module load BCFtools/1.15-GCC-11.2.0

cat $RESULTS/plot*.txt | awk -v thres=$Threshold  '$4 < 0.05/thres'|awk '{print $1"\t"$3}'|while read line
	do Chr=$(echo $line|awk '{print $1}'); SNP=$(echo $line|awk '{print $2}')
	   paste <(cat $PHENOTYPES) <(bcftools view --regions $Chr:"$SNP" $VCF| bcftools query -f '%POS  [ %GT]\n'] |perl -pe 's/ /\n/g'|awk 'NR > 3') <(cat $SPECIES |awk -v CHR=$Chr -v POS=$SNP '{print $2"\t"CHR"\t"POS}') |grep -v -E "\-9"|grep -P "\t2\t"
	   paste <(cat $PHENOTYPES) <(bcftools view --regions $Chr:"$SNP" $VCF| bcftools query -f '%POS  [ %GT]\n'] |perl -pe 's/ /\n/g'|awk 'NR > 3') <(cat $SPECIES |awk -v CHR=$Chr -v POS=$SNP '{print $2"\t"CHR"\t"POS}') |grep -v -E "\-9"|grep -P "\t3\t"
     paste <(cat $PHENOTYPES) <(bcftools view --regions $Chr:"$SNP" $VCF| bcftools query -f '%POS  [ %GT]\n'] |perl -pe 's/ /\n/g'|awk 'NR > 3') <(cat $SPECIES |awk -v CHR=$Chr -v POS=$SNP '{print $2"\t"CHR"\t"POS}') |grep -v -E "\-9"|grep -P "\t2.5\t"
	done 

