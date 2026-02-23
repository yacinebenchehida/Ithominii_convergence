
# Transpolymorphism Sliding Window Pipeline

This pipeline identifies transpolymorphic sites across multiple species from a VCF file, summarizes their distribution in sliding windows, and produces genome-wide diagnostic plots.

The workflow is controlled by:

- SLURM wrapper script (job submission)
- `transP_v3.py` (core analysis and plotting)

---

## Required Software

- Python ≥ 3.8
- pysam
- matplotlib
- BCFtools (for chromosome extraction)

---

## Input Files

### 1. VCF File

- bgzipped VCF (.vcf.gz)
- Must contain **invariant** and variant sites
- Must include all samples listed in the population file

### 2. Population File

Tab-separated file:

sampleID    species

Used to group samples by species and compute species-level allele counts.

---

## Script Usage

```
python3 transP_v3.py   -v <VCF>   -p <POP_FILE>   -c <CHROMOSOME>   -o <OUTPUT_PREFIX>   -w <WINDOW_SIZE_BP>   -d <WORKING_DIRECTORY>   --highlight_start <START_BP>   --highlight_end <END_BP>
```

---

## Pipeline Steps

### 1. Chromosome Extraction

Uses:

```
bcftools view -r <chromosome>
```

Creates:
```
extracted_chromosome.vcf.gz
```

Stored inside the working directory.

---

### 2. Species-Level Allele Counting

For each VCF record:

- Samples are grouped by species
- Alleles are counted per species
- Missing genotypes are ignored

Produces a dictionary:

species → {allele: count}

---

### 3. Transpolymorphic Site Definition

A site is considered:

**Polymorphic in a species** if ≥ 2 alleles are present.

If number of species > 3:
- At least 4 species must be polymorphic

If number of species ≤ 3:
- All species must be polymorphic

Additional allele count filter:
- Every allele within every species must have count ≥ 6

Only sites passing all criteria are retained as transpolymorphic.

The script tracks:

- Total sites
- All polymorphic sites
- Transpolymorphic sites

---

### 4. Sliding Window Analysis

Window size is defined in base pairs (`-w`).

For each window:

- Count transpolymorphic sites
- Count available sites
- Count polymorphic sites

Compute:

- Ratio to available sites
- Ratio to polymorphic sites

Outputs:

```
window_counts.txt
```

---

### 5. Plot Generation

Three PDF plots are generated:

1. `transpolymorphic_sites_plot.pdf`
   - Number of transpolymorphic sites per window

2. `available_sites_ratio_plot.pdf`
   - Ratio of transpolymorphic sites to total available sites

3. `polymorphic_sites_ratio_plot.pdf`
   - Ratio of transpolymorphic sites to polymorphic sites

The specified genomic interval (`--highlight_start`, `--highlight_end`) is shaded in red.

---

## Output Structure

Inside the working directory:

- extracted_chromosome.vcf.gz
- <OUTPUT_PREFIX> (transpolymorphic site list)
- window_counts.txt
- transpolymorphic_sites_plot.pdf
- available_sites_ratio_plot.pdf
- polymorphic_sites_ratio_plot.pdf

---

## Notes

- No genotype-level quality filtering is performed.
- Allele count threshold (≥6 per allele per species) is hard-coded.
- Windowing is coordinate-based (fixed bp size).
- Requires sufficient per-species sample depth to pass allele count threshold.
