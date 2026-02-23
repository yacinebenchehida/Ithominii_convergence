# f4 / f3 Statistics Pipeline (ADMIXTOOLS)

This pipeline computes **f4 statistics** (primary use) and optionally **f3 statistics** using the R package `admixtools`, starting from a multi-sample VCF.

The workflow:

1. Subset a VCF to specific populations
2. Filter SNPs (biallelic, low missingness, GQ >= 10)
3. Convert to PLINK format
4. Run `admixtools::f4()` (primary use)
5. Write tab-delimited results

The main execution script is:

- `f4_f3.sh` (submitted via `sbatch`)

The master script simply loops over multiple population combinations and submits jobs.

---

## Files and Roles

### `master` (loop script)
Defines:
- Input VCF
- Output directory
- Population file
- Lists of P1–P4 combinations

Submits jobs:

```bash
sbatch ./f4_f3.sh   -v <vcf>   -o <output_dir>   -p <pop_file>   --p1 <pop1>   --p2 <pop2>   --p3 <pop3>   --p4 <pop4>   --stat f4
```

---

### `f4_f3.sh`
Main SLURM job script. Performs:

1. Argument parsing
2. Population selection
3. VCF subsetting and filtering
4. PLINK conversion
5. Execution of `f4_f3.R`

---

### `f4_f3.R`
Calls:

- `admixtools::f4()` for f4
- `admixtools::qp3pop()` for f3

Writes results as TSV files.

---

# 1) Inputs

## 1.1 VCF (`-v`)

- Must be **bgzipped and indexed**
- Multi-sample VCF
- All populations listed in the population file must be present

Filtering applied during subsetting:

- `F_MISSING < 0.2`
- `FORMAT/GQ >= 10`
- Biallelic SNPs only (`-m2 -M2 -v snps`)

---

## 1.2 Population File (`-p`)

Two-column tab-delimited file:

```
population_name    sample_ID
```

Column 1: population label  
Column 2: sample ID (must match VCF)

Selection is done using:

```
awk '$1 ~ pattern'
```

This means population labels are matched using regex pattern matching.

---

# 2) Parameters

## Required

- `-v, --vcf`      Input VCF (bgzipped)
- `-o, --outdir`   Output directory
- `-p, --pop`      Population file
- `--p1`           Population 1
- `--p2`           Population 2
- `--p3`           Population 3
- `--stat`         `f4` or `f3`

## Required for f4

- `--p4`           Population 4

Note: Although f3 is implemented, this pipeline is primarily designed for **f4 analyses**.

---

# 3) Pipeline Logic

## Step 1: Create phenotype and subspecies files

Two files are generated in:

```
<OUTDIR>/<PREFIX>/
```

- `<PREFIX>_phenotype_file.txt`
- `<PREFIX>_subspecies.txt`

`PREFIX`:

- For f4: `P1_P2_P3_P4`
- For f3: `P1_P2_P3`

---

## Step 2: Subset and Filter VCF

Using:

```bash
bcftools view   -S <sample_list>   -i 'F_MISSING<0.2 && FORMAT/GQ>=10'   -m2 -M2 -v snps
```

Output:

```
<PREFIX>.vcf.gz
```

---

## Step 3: Convert to PLINK

```
plink --vcf <vcf> --make-bed
```

Produces:

- `Inputs.bed`
- `Inputs.bim`
- `Inputs.fam`

Contigs are converted to numeric IDs to ensure compatibility with ADMIXTOOLS.

---

## Step 4: Compute f4 (Primary Use)

In `f4_f3.R`:

```r
out <- f4(data_plink, pop1, pop2, pop3, pop4)
```

This computes:

```
f4(pop1, pop2; pop3, pop4)
```

Using the `admixtools` R implementation.

The output contains:

- f4 estimate
- standard error
- Z-score

Block jackknife is handled internally by `admixtools`.

Results are written as:

```
pop1_pop2_pop3_pop4_f4.txt
```

---

# 4) Interpretation of f4

f4 tests the symmetry of allele frequency correlations:

```
f4(A, B; C, D)
```

Interpretation:

- f4 ≠ 0 → allele frequency asymmetry
- |Z| > 3 typically considered significant
- Sign indicates direction of correlation

Positive or negative values reflect excess allele sharing between specific pairs.

---

# 5) Optional f3

If `--stat f3`:

```r
qp3pop(data_plink, pop1, pop2, pop3)
```

Computes:

```
f3(pop1; pop2, pop3)
```

Negative f3 with significant Z-score indicates admixture in pop1.

Output file:

```
pop1_pop2_pop3_f3.txt
```

---

# 6) Output Structure

For each test:

```
<OUTDIR>/<PREFIX>/
    <PREFIX>_phenotype_file.txt
    <PREFIX>_subspecies.txt
    <PREFIX>.vcf.gz
    Inputs.bed
    Inputs.bim
    Inputs.fam
    pop1_pop2_pop3_pop4_f4.txt   (or f3 equivalent)
```

---

# 7) Software Requirements

Modules loaded in SLURM script:

- BCFtools
- PLINK 1.9
- R 4.2+
- R package: admixtools
- R package: tidyverse

---

# 8) Notes

- Designed primarily for **f4 statistics**
- f3 implemented for convenience only
- Uses genotype quality filtering (GQ >= 10)
- Excludes sites with >20% missing data
- Uses only biallelic SNPs
- No LD pruning is performed
- No MAF filtering is applied

---
