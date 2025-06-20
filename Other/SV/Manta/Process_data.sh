module load BCFtools/1.15.1-GCC-11.3.0

# Extract SV genotypes (separating ssp and hicetas)
for i in {1..16}; do bcftools view --regions SUPER_18:15200397-15263591 diploidSV.vcf.gz|bcftools  query -f '%POS\t%REF\t%ALT\t%INFO/SVLEN  [ %GT]\n']|awk -v SV=$i 'NR==SV'|cut -f 1-3; echo ssp; cat list_ssp.txt|while read line; do paste <(bcftools query -l diploidSV.vcf.gz) <(bcftools view --regions SUPER_18:15200397-15263591 diploidSV.vcf.gz|bcftools  query -f '%POS\t%ALT\t%INFO/SVLEN  [ %GT]\n']|awk -v SV=$i 'NR==SV'|perl -pe 's/(\s)+/\n/g'|awk 'NR > 3')|grep $line; done; echo hicetas;cat list_hicetas.txt|while read line; do paste <(bcftools query -l diploidSV.vcf.gz) <(bcftools view --regions SUPER_18:15200397-15263591 diploidSV.vcf.gz|bcftools  query -f '%POS\t%ALT\t%INFO/SVLEN  [ %GT]\n']|awk -v SV=$i 'NR==SV' |perl -pe 's/(\s)+/\n/g'|awk 'NR > 3')|grep $line; done; done

# Look specifically to the genotype at deletions (separating ssp and hicetas)
for i in {1..11}; do bcftools view diploidSV.vcf.gz |bcftools  query -f '%POS\t%REF\t%ALT\t%INFO/SVLEN  [ %GT]\n'|grep -P "<DEL>"|awk -v SV=$i 'NR==SV'|cut -f 1-3; echo ssp; cat list_ssp.txt|while read line; do paste <(bcftools query -l diploidSV.vcf.gz) <(bcftools view diploidSV.vcf.gz|bcftools  query -f '%POS\t%ALT\t%INFO/SVLEN  [ %GT]\n']|awk -v SV=$i 'NR==SV'|perl -pe 's/(\s)+/\n/g'|awk 'NR > 3')|grep $line; done;\
 echo hicetas; cat list_hicetas.txt|while read line; do paste <(bcftools query -l diploidSV.vcf.gz) <(bcftools view diploidSV.vcf.gz|bcftools  query -f '%POS\t%ALT\t%INFO/SVLEN  [ %GT]\n']|awk -v SV=$i 'NR==SV'|perl -pe 's/(\s)+/\n/g'|awk 'NR > 3')|grep $line; done;done

