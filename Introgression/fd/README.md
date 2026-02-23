
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

---

## Plots

### 1. All statistics plot

Plot showing:

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

---
