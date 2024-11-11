module load  BCFtools/1.19-GCC-13.2.0
module load VCFtools/0.1.16-GCC-11.2.0
module load Pysam/0.22.0-GCC-12.3.0

VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/multisp/Melinaea_menophilus/multisp_Melinaea_menophilus_genotypeGVCF.intervals_55.filters.DP4_GQ5_QUAL5.invariants.vcf.gz"
CHR="SUPER_18"
START="15090001"
END="16000000"
POP="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/MuteBass/Inputs/yellow_band_species.txt"

# Ingroup
bcftools view --regions $CHR:$START-$END -S <(cat $POP|awk '{print $1}') $VCF | bcftools view -U -i 'TYPE=="snp"' -m2 -M2 -v snps -O v|bcftools view -i 'F_MISSING < 0.5'  > test.vcf
bgzip test.vcf
tabix test.vcf.gz

# outgroup
bcftools view --regions $CHR:$START-$END -T <(bcftools query -f '%CHROM\t%POS\n' test.vcf.gz) -S <(cat $POP|awk '$2=="outgroup" {print $1}') $VCF > outgroup.vcf
bgzip outgroup.vcf
tabix outgroup.vcf.gz

# Allele frequencies outgroup
indv_command=$(for i in $(bcftools query -l  outgroup.vcf.gz); do echo -e "--indv $i"; done |perl -pe 's/\n/ /g')
vcftools --gzvcf outgroup.vcf.gz --freq $indv_command --out outgroups


awk '{
    split($5, ref, ":"); 
    split($6, alt, ":"); 
    if (ref[2] == "-nan") { 
        print $2, "nan"; 
    } else if (ref[2] > 0.5) { 
        print $2, ref[1]; 
    } else if (ref[2] == 0.5) { 
        print $2, "nan"; 
    } else { 
	print $2, alt[1]; 
    }
}' <(grep -v CHROM outgroups.frq) > ancestral.alleles

python3 select_random.py $POP

for i in $(cat $POP|awk '{print $2}'|sort -u)
do
	bcftools view --regions $CHR:$START-$END  -T <(bcftools query -f '%CHROM\t%POS\n' test.vcf.gz) -S <(cat "$i".txt|awk '{print $1}') $VCF > "$i".vcf
	bgzip "$i".vcf
	tabix "$i".vcf.gz
	vcftools --gzvcf "$i".vcf.gz --freq --out "$i".frequencies
	grep -f <(cat ancestral.alleles|grep nan |awk '{ print $1}') "$i".frequencies.frq >  "$i".frq
	bcftools query -f '%POS\t%REF\t%ALT\t[%GT\t]\n' "$i".vcf.gz > "$i".geno
	rm "$i".frequencies.frq
done

SPECIES=$(cat $POP|awk '{print $2}'|grep -v outgroup|sort -u|perl -pe 's/\n/,/g')
echo $SPECIES


module load R/4.2.1-foss-2022a
module load Julia/1.9.3-linux-x86_64
Rscript ./data_table_version.r tarapotensis,mothone,menophilus,isocomma,marsaeus ./

#python3 ./input_mutebass.py -a ancestral -v test.vcf.gz -s /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/MuteBass/Inputs/yellow_band_species.txt
