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
### Generate GVCF
- The folder 2_gvc contains the scripts that were used to obtain the GVCFs using GATK HaplotypeCaller. For the sake of speed, the script was applied to each genome and each genome interval.

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
  


## Example usage
### Run the full pipeline
``` bash
library(SAMPLE)
data("coral_symbionts")
set.seed(812)
SAMPLE(input = coral_symbionts, output_N = "Example", replicates = 50, stability_thresh = 2, sucess_points = 10, diff = 1)
```
