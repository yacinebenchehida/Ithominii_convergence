# TWISST launcher (launch_twisst)

This folder contains a pipeline to:

1. Subset a multi-sample VCF to a region and a set of taxa
2. Split the region into **windows containing a constant number of variant SNPs** (while retaining invariant sites)
3. Infer a phylogeny for each window (NJ or ML)
4. Combine all window trees into a single tree file
5. Run **TWISST** topology weighting
6. Plot topology weights and optionally run a block permutation significance test if there is an enrichment of the introgression topology within the peak

---

## Files and roles

- `launch_twisst`  
  Main launcher: subsets VCF, creates phenotype file, splits into window VCFs, submits per-window phylogeny jobs, then submits TWISST as a dependent job.

- `splice_VCF_constant_SNPs.py`  
  Splits a bgzipped VCF into multiple VCF parts, each containing **SNPS** variant sites (ALT != "."). Invariant sites are retained. Logs the genomic range (start/end position) for each part.

- `Phylogenies_job_submitter.py`  
  Creates and submits SLURM job scripts to convert each window VCF to PHYLIP and infer a phylogeny for each window using:
  - `NJ` (R/ape via `NJ_tree.R`) or
  - `ML` (RAxML-NG)

- `run_twisst.sh`  
  Combines all per-window trees, removes windows with missing trees, runs TWISST, calls `twisst.R` for plotting and permutations, and cleans temporary files.

- `twisst.R`  
  Produces TWISST plots and, when peak boundaries are provided, runs a block permutation test and writes `significance.txt`.

---

## 0) Running the pipeline

Submit with `sbatch`:

```bash
sbatch ./launch_twisst \
  -v <vcf> \
  -c <chromosome> \
  -s <start_position> \
  -e <end_position> \
  --snps <snps_per_window> \
  -f <pop_file> \
  -t <list_taxa> \
  -m <method> \
  -o <path_output> \
  -n <name> \
  --introg <taxon1,taxon2> \
  --ps <peak_start> \
  --pe <peak_end>
```

### Example

```bash
sbatch ./launch_twisst \
  -v /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/multisp/Melinaea_marsaeus/multisp_Melinaea_marsaeus_genotypeGVCF.intervals_8.filters.DP4_GQ5_QUAL5.invariants.vcf.gz \
  -c SUPER_2 \
  -s 25565242 \
  -e 26521063 \
  --snps 50 \
  -f /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/Relate/Inputs/hindwing_black.txt \
  -t mothone,satevis,marsaeus_rileyi,marsaeus_phasiana \
  -m NJ \
  -o ./ \
  -n test_v2 \
  --introg marsaeus_rileyi,mothone \
  --ps 25965242 \
  --pe 26021063
```

---

## 1) Inputs

### 1.1 VCF (`-v`)

- Input VCF file name (typically `.vcf.gz`)
- Must contain the samples referenced in the population file
- Must contain the requested region (`-c`, `-s`, `-e`)

### 1.2 Population file (`-f`)

Tab-delimited file used to select samples matching the taxa in `-t`.  
The launcher writes a phenotype file containing 2 columns:

- Column 1: sample ID
- Column 2: group (used by TWISST via `--groupsFile`)

`launch_twisst` builds this file using:

- `grep -P "\\s<i>\\s*(\\w)*\\s*$"` for each taxon label `i` in `-t`
- then `awk '{print $1"\\t"$2}'`

This means the population file must contain the taxon label as a **whole field** near the end of each matching line, and the first two columns must be sample ID and group.

---

## 2) Parameters

### Mandatory parameters

- `-v <vcf>`  
  Input VCF file name

- `-c <chromosome>`  
  Contig / scaffold / chromosome name

- `-s <start_position>`  
  Start position in the VCF region

- `-e <end_position>`  
  End position in the VCF region

- `--snps <snps_per_window>`  
  Number of **variant SNPs** per window. Each output window contains this many sites where ALT != "." (invariant sites are still written, so window lengths in bp vary).

- `-f <pop_file>`  
  Population file

- `-t <list_taxa>`  
  Comma-separated taxa labels to include (used to select samples)

- `-m <method>`  
  Phylogeny method:
  - `NJ` (R/ape)
  - `ML` (RAxML-NG)

- `-o <path_output>`  
  Output directory

- `-n <name>`  
  Prefix for output files and output subdirectory (`<path_output>/<name>`)

- `--introg <taxon1,taxon2>`  
  Two taxa names used to define the focal introgression topology in `twisst.R`  
  (required by the launcher script)

### Optional parameters

- `--ps <peak_start>` and `--pe <peak_end>`  
  Peak boundaries. When provided, `twisst.R`:
  - zooms the plots to a 500 kb window around the peak midpoint
  - draws the peak as a rectangle
  - runs permutation tests and writes `significance.txt`

---

## 3) Pipeline steps (internal logic)

### Step 1: Create working directory and phenotype file

The launcher creates:

- `<output_path>/<name>/`
- `<name>_phenotype_file.txt` inside that directory

This file is used as:
- the sample list for `bcftools view -S`
- the TWISST `--groupsFile` (via `run_twisst.sh`)

### Step 2: Subset the VCF to region and samples

`bcftools view` is used to subset by:
- region: `CHR:START-END`
- samples: from `<name>_phenotype_file.txt`

Output:
- `<name>.vcf.gz` (bgzipped and indexed)

### Step 3: Split the region into constant-SNP windows

`splice_VCF_constant_SNPs.py` creates:

- `<name>_part_1.vcf`, `<name>_part_2.vcf`, ...
- `<name>.txt` log file with 3 columns:

```text
part_id    start_position    end_position
```

Important detail: each part contains **SNPS** variant sites (ALT != "."); invariant sites are included in the output VCF parts.

### Step 4: Infer a phylogeny per window (SLURM batch jobs)

`Phylogenies_job_submitter.py`:

- groups VCF parts into batches of 1000 per SLURM script (`job_1.sh`, `job_2.sh`, ...)
- for each window VCF part:
  - bgzip + tabix
  - convert VCF to PHYLIP with `vcf2phylip.py -r -m 0`
  - infer a tree using:
    - NJ: `Rscript NJ_tree.R ...` and retain only `*.newick`
    - ML: `raxml-ng-mpi --search ...` and retain only `*.bestTree`
  - deletes the window VCF and its index after processing

Tree outputs are written into per-window directories under `<output_path>/<name>/`, named by the VCF part prefix.

### Step 5: Run TWISST as a dependent job

The launcher:

1. queries running SLURM jobs matching a derived `job_name`
2. submits `run_twisst.sh` with:

```bash
sbatch --dependency=aftercorr:<jobids> --job-name="twi<job_name>" ./run_twisst.sh ...
```

`aftercorr` ensures TWISST runs after the phylogeny jobs complete.

---

## 4) TWISST execution and outputs (`run_twisst.sh`)

### 4.1 Combine per-window trees

Inside `<output_path>/<name>/`:

- concatenates all available per-window tree files into:
  - `topologies.txt`

Tracks missing windows in:
- `missing_parts.txt`

### 4.2 Build positions file

From `<name>.txt` (the window range log):

- if there are missing windows, those rows are removed
- writes:

```text
positions.txt
```

with 2 columns:

```text
start    end
```

### 4.3 Run TWISST

Creates `-g` options from the phenotype file (unique groups in column 2), then runs:

- `Phylogeny.topos` (written by TWISST via `--outputTopos`)
- `Phylogeny.weights.tsv.gz` (weights per window)

### 4.4 Plot and permutations

Runs:

```bash
Rscript twisst.R <start> <end> 0.012 <taxon1> <taxon2> <peak_start> <peak_end>
```

Outputs produced by `twisst.R` include:

- `Contribution.pdf`
- `Principal_topology.pdf`
- `twisst_barplot_no_smoothing.pdf`
- `twisst_barplot_with_smoothing.pdf`

If peak boundaries are provided, also:

- `significance.txt`

### 4.5 Cleanup

`run_twisst.sh` removes:

- all window directories: `<name>_part_*`
- `*sh` in the output directory

---

## 5) Permutation test (implemented in `twisst.R`)

When peak boundaries are provided (peak mode), `twisst.R` performs a block permutation test on the focal introgression topology.

### 5.1 Focal topology definition

- Reads `Phylogeny.topos`
- Finds the topology index where the two taxa (`taxon1`, `taxon2`) are siblings
- Extracts that topology weight from `Phylogeny.weights.tsv.gz` as `introgression_topo`

A window is considered high support if:

```text
introgression_topo >= 0.95
```

Observed statistic:

- number of windows within `[peak_start, peak_end]` with `introgression_topo >= 0.95`

### 5.2 Permutation schemes

Parameters:

- `n_permutations = 50000`
- `block_size = 100000` bp

Two schemes are run:

1. **Inter-block permutation (between blocks only)**  
   Blocks are permuted; within-block order is preserved.

2. **Inter + intra-block permutation (between blocks + within-block shuffling)**  
   Blocks are permuted; within each block, `introgression_topo` values are shuffled.

### 5.3 Reported outputs

`significance.txt` with columns:

- `zscore_inter` and `pvalue_inter`
- `zscore_inter_intra` and `pvalue_inter_intra`

If the observed statistic is 0:

- z-scores and p-values are recorded as:
  - `No intro Topo > 95% in peak`

---

## 6) Expected directory structure (after completion)

```text
<output_path>/<name>/
  <name>_phenotype_file.txt
  <name>.vcf.gz
  <name>.vcf.gz.tbi
  <name>.txt
  topologies.txt
  missing_parts.txt
  positions.txt
  Phylogeny.topos
  Phylogeny.weights.tsv.gz
  Contribution.pdf
  Principal_topology.pdf
  twisst_barplot_no_smoothing.pdf
  twisst_barplot_with_smoothing.pdf
  significance.txt              (only if peak mode)
```

Note: per-window directories (`<name>_part_*`) are removed at the end by `run_twisst.sh`.

---
