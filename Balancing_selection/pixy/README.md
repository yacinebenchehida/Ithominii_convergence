
# pixy Nucleotide Diversity Pipeline

This pipeline computes nucleotide diversity (π) in sliding windows using `pixy`, applies region masking, and produces a genome-wide plot highlighting a peak interval.

---

## Required Software

- BCFtools
- pixy
- R (ggplot2)

---

## Input Files

### 1. VCF
- bgzipped and indexed (.vcf.gz)
- Contains all samples of interest

### 2. Population file
Tab-separated, two columns:

sampleID    population

Samples with population value `-9` are excluded automatically.

### 3. Mask file
Tab-separated file:

chromosome    start    end

Regions matching the target chromosome and interval are excluded before running pixy.

---

## What the Script Does

1. Creates output directory.
2. Filters population file (removes `-9` entries).
3. Extracts target genomic interval from VCF.
4. Removes masked regions.
5. Runs:

   pixy --stats pi  
   --window_size 10000  
   --chromosomes <CHR>

6. Reformats pixy output for plotting.
7. Generates a PDF plot with the peak interval highlighted (using plot_nucldiv.R)

---

## Outputs

Inside the specified results directory:

- `<prefix>_pi.txt`  
  Raw pixy output

- `<prefix>_pi2plot.txt`  
  Midpoint + π value per window

- `<prefix>_nucleotide_diversity.pdf`  
  Genome-wide π plot with peak region highlighted
