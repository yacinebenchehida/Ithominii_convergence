# Mapping, SNPs calling and phasing: 

This folder contains all the scripts used to perform all the step from the initial mapping of the raw reads to the phasing of the final VCF.

## 0) Preparing the reference genomes for mapping and snp calling

Before use, fasta reference genome were:
- Indexed for bwa:
``` bash
bwa index $ref_genome
```
- Indexed for samtools ("faidxed")
``` bash
samtools faidx $ref_genome
```
- Created a GATK dictionary 
``` bash
gatk CreateSequenceDictionary -R $ref_genome
```
- Split into 60 intervals to speed up downstream analyses
``` bash
gatk SplitIntervals -R $ref_genome --scatter-count 60 -O 
```

## 1) Mapping

- The folder 1_mapping contains the scripts that were used to map the reads to the reference genome using bwa mem, transform the sam into a bam file and then sort the reads in the bam file:

``` bash
bwa mem -t 12 $ref_genome R1.fastq.gz R2.fastq.gz | samtools view -bSh |samtools sort -o sorted.bam
```

- Then the Picard was used to add reads groups and mark duplicates

``` bash
# Mark read groups
java -jar -Xmx50g $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
 I=$path_results/0_sorted_bam/$1/"$1"_sorted.bam \
 O=$path_results/0_sorted_bam/$1/"$1"_sorted_RG.bam \
 RGID=$RGID \
 RGLB=$RGLB \
 RGPL=$RGPL \
 RGPU=$RGPU \
 RGSM=$RGSM \
 SORT_ORDER=coordinate \
 CREATE_INDEX=TRUE

# Mark duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
VALIDATION_STRINGENCY=SILENT \
I=sorted_RG.bam \
O=sorted_dedup.bam \
M=marked_dedup_metrics.txt
```

- Finally Qualimap was used to summarize and inspect the quality of the mapping
``` bash
qualimap bamqc --java-mem-size=20G \
-bam sorted_dedup.bam \
-c \
-gd HUMAN \
-nt 12 \
-outformat HTML 
```

## 2) SNPs calling
The SNP calling was performed in several steps:

### Generate GVCF
The folder 2_gvcf contains the scripts that were used to obtain the GVCFs using GATK HaplotypeCaller. For the sake of speed, the script was applied to each genome and each genome interval.

``` bash
gatk --java-options "-Xmx4g" HaplotypeCaller \
-R $ref_genome \
-I sorted_dedup.bam \
-O ${SLURM_ARRAY_TASK_ID}.g.vcf.gz \
--intervals $intervals \
--output-mode EMIT_ALL_CONFIDENT_SITES \
-ERC GVCF \
--dont-use-soft-clipped-bases 
```
  
### Combine GVCFs
The folder 3_Combined_gvcf contains the scripts used to merge GVCFs of all samples of one species. We used the CombineGVCFs option of GATK for that end. For the sake of speed, the script was applied to each species and each genome interval.

``` bash
gatk --java-options -Xmx15g CombineGVCFs \
-R $ref_genome \
$variants \
-O intervals_merged_${SLURM_ARRAY_TASK_ID} \
-intervals $intervals
```

### Genotype the GVCFs
The folder 4_Genotype_gvcfs contains the scripts necessary to genotype the GVCFs and generate raw VCFs. It uses GenotypeGVCFs option of GATK. For the sake of speed, the script was applied to each species and each genome interval.

``` bash
gatk --java-options -Xmx15g GenotypeGVCFs \
-R $ref_genome \
-V intervals_merged_${SLURM_ARRAY_TASK_ID} \
-O genotypeGVCF.intervals_${SLURM_ARRAY_TASK_ID}.vcf.gz
-intervals $intervals
```
### Filter VCFs
The VCFs were further filtered using BCFtools. More precisely BCFtools was used to keep bi-allelic SNPs with a variant quality score of at least 10, a genotype quality of at least 10, a depth of coverage of at least 5 and to exclude any SNPs with more than 20% of missing data.

``` bash
bcftools filter -e 'FORMAT/DP < 1 |FORMAT/GQ < 5 |QUAL <= 5' --set-GTs . genotypeGVCF.intervals_${SLURM_ARRAY_TASK_ID}.vcf.gz -O u | bcftools view -U -i 'TYPE=="snp"' -m2 -M2 -v snps -O v| bcftools view -i 'F_MISSING < 0.2'> genotypeGVCF.intervals_${SLURM_ARRAY_TASK_ID}.filters.snps.vcf
```
### Combining all VCF intervals in a single VCF
We combined all 60 intervals to generate a single VCF for each species. Each VCF was also zipped (bgzip) and tabulated (tabix). 

``` bash
cat *intervals_0.* > snps.vcf
for j in {1..59}; do cat *intervals_"$j".*|grep -v '#' >> snps.vcf; done
bgzip -@ 12 snps.vcf && tabix -p vcf snps.vcf.gz
```

## 3) Phasing and missing data imputation
The phasing and the data imputation were done using ShapeIT with default option. Shape it was applied to each scaffold/chromosome separately. Each phased VCF was also zipped (bgzip) and tabulated (tabix). 

``` bash
SCAFFOLD=$(sed -n "${SLURM_ARRAY_TASK_ID}"p list_of_scaffolds.txt)
$SHAPEIT --input snps.vcf.gz --region $SCAFFOLD --output phased.snps.vcf
bgzip phased.snps.vcf && tabix -p vcf phased.snps.vcf.gz
```


