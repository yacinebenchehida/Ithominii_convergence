
# fd Sliding Window Pipeline

## Overview

This pipeline computes **D-statistic, fd, and fdM** in sliding windows using the [genomics_general](https://github.com/simonhmartin/genomics_general) toolkit (ABBABABAwindows.py).  
It performs:

1. VCF subsetting to selected taxa
2. SNP filtering
3. Conversion to `.geno.gz` format
4. Sliding window ABBA–BABA statistics
5. Filtering of low-quality windows
6. Visualization of D, fd, fdM and ABBA/BABA counts

---

## Dependencies

The following software must be available:

- HTSlib
- BCFtools
- R (≥ 4.2)
- Python3
- `genomics_general` toolkit  
  - `parseVCF.py`
  - `ABBABABAwindows.py`

The script assumes the path to [genomics_general](https://github.com/simonhmartin/genomics_general) is correctly defined inside the `fd` script.

---

## Input Files

### 1. VCF file

- Must be **bgzipped and indexed**
- Must contain all samples referenced in the population file

### 2. Population file (`pop.txt`)

Tab-separated file with at least:

Column 1: Sample ID  
Column 2: Population name  

The population names must match the values provided to:

- `-p1`
- `-p2`
- `-p3`
- `-po`

---

## fd Script Usage

```
./fd  -v <vcf>       -c <chromosome>       -s <start_position>       -e <end_position>       -f <pop_file>       -p1 <taxa1>       -p2 <taxa2>       -p3 <taxa3>       -po <outgroup>       -w <window_size>       -o <output_path>       -n <prefix>
```

### Parameters

- `-v`   Input VCF (bgzipped)
- `-c`   Chromosome/contig name
- `-s`   Start coordinate
- `-e`   End coordinate
- `-f`   Population file
- `-p1`  Population 1
- `-p2`  Population 2
- `-p3`  Population 3
- `-po`  Outgroup population
- `-w`   Window size (in number of SNPs)
- `-o`   Output directory
- `-n`   Prefix for output files

---

## Pipeline Steps

### 1. Population Subsetting

Creates a phenotype file containing samples for:

P1, P2, P3, and Outgroup.


---

### 2. Convert VCF to GENO

Using:

```
parseVCF.py
```

Generates:

```
<PREFIX>.geno.gz
```

---

### 3. Sliding Window ABBA–BABA

Using:

```
ABBABABAwindows.py
```

Options used:

- `--windType coordinate`
- `-w <window_size>`
- `-m 1` (minimum 1 SNP per window)
- `-f phased`
- `--writeFailedWindows`

Statistics computed per window:

- ABBA
- BABA
- D
- fd
- fdM
- sitesUsed

Output:

```
output.csv
```

Converted to tab-delimited:

```
Results_<CHR>_<START>_<END>_<P1>_<P2>_<P3>_<O>.txt
```

Windows with `nan` are removed.

---

## Data Filtering (Statistics.R)

The following filters are applied before plotting:

- D < 0 → set to 0
- fd < 0 → set to 0
- fd > 1 → set to 0
- fdM < 0 → set to NA
- Remove lowest 2% of windows by `sitesUsed`

This removes poorly supported windows.

---

## Output Files

Inside:

```
<OUTPUT_PATH>/<PREFIX>/
```

Generated files include:

- `<PREFIX>_phenotype_file.txt`
- `<PREFIX>.vcf.gz`
- `<PREFIX>.geno.gz`
- `output.csv`
- `Results_*.txt`
- `<PREFIX>_All_stats_plot.pdf`
- `<PREFIX>_ABBA_BABA_stats_plot.pdf`

---

## Plots

### 1. All statistics plot

Faceted plot showing:

- ABBA
- BABA
- D
- fd
- fdM

Across genomic coordinates.

### 2. ABBA/BABA plot

Shows only:

- ABBA
- BABA

With title indicating P1, P2, P3, and Outgroup.

---

## Interpretation Notes

- Positive D indicates excess ABBA over BABA.
- fd and fdM quantify introgression proportion.
- Windows are defined by fixed SNP count (`-w`), not fixed physical size.
- fd values are constrained to [0,1] after filtering.
- Low coverage windows (lowest 2%) are removed.

---

## Master Script

The master script submits multiple combinations of P1, P2, P3 using:

```
sbatch ./fd ...
```

Each run produces an independent results directory.

---

## Assumptions

- VCF is correctly filtered upstream.
- Population names in pop.txt are exact matches.
- SNPs are bi-allelic.
- Data are phased (even though phasing step is currently commented).

---

## Summary

This pipeline computes sliding-window ABBA–BABA statistics and produces:

- Genome-wide D
- fd
- fdM
- Diagnostic ABBA/BABA counts
- Publication-ready PDF plots

It is designed for high-throughput testing of many population combinations on HPC clusters.
