# RELATE pipeline (Relate_launcher_V4)

This folder contains all scripts used to run the **RELATE** analyses, including phasing (SHAPEIT4), ancestral state inference, SNP specific tree extraction, and TWISST topology weighting.

## Contents

- `Relate_launcher_V4` (main launcher, submit with `sbatch`)
- Helper scripts (must be in the same folder as the launcher):
  - `createuniformrecmap.r`
  - `Find_Alternative_ancestral.r`
  - `rename_tree.py`
  - `tree_extract.py`
  - `plot_tree.R`
  - `position_4_twisst.py`
  - `twisst.R`
  - `twisst_permutations.R`

---

## 0) Running the whole pipeline

Submit the pipeline with `sbatch`:

```bash
sbatch ./Relate_launcher_V4 \
  -v <vcf> \
  -c <chromosome> \
  -s <start_position> \
  -e <end_position> \
  --snps <list_snps> \
  -f <pop_file> \
  -t <list_taxa> \
  -r <root_outgroup> \
  --species <species_list> \
  -o <path_output> \
  -n <name> \
  --introg <two_species_names> \
  --ps <peak_start> \
  --pe <peak_end>
```

### Example

```bash
sbatch ./Relate_launcher_V4 \
  -v /path/to/My_VCF.vcf.gz \
  -c SUPER_4 \
  -s 1 \
  -e 1843557 \
  -f ../Inputs/yellow_band_info.txt \
  -t mothone_mothone,mothone_messenina,isocomma,simulator \
  --snps 1,1390852,1390869,1391864,1391929,1391932,1392037,1392064,1392096,1392107,1392262,1392265,1392271,1392274,1392382,1392390,1392438,1392491,1392714,1395239,1398720,1843557 \
  -r outgroup \
  --species tarapotensis,satevis,lilis,isocomma,idae,marsaeus,menophilus,mothone,flavo,mneme,maeonis \
  -o ../Results/Melinaea_mothone_ref/introgression_yellow_band \
  -n mothone_messenina_isocomma_simulator \
  --introg mothone_messenina,simulator \
  --ps 1385004 \
  --pe 1398720
```

---

## 1) Inputs

### 1.1 VCF (`-v`)

- Must be bgzipped and indexed (`.vcf.gz` with a `.tbi`)
- Must contain all samples referenced in the population file
- Must contain the requested region (`-c`, `-s`, `-e`)

### 1.2 Population file (`-f`)

Tab delimited file containing at least two columns:

- Column 1: sample ID (must match the VCF sample IDs)
- One of the later columns must contain the taxon labels used by `-t`, `-r`, and `--species`

Example:

```text
sample1	mothone_mothone
sample2	mothone_mothone
sample3	mothone_messenina
sample4	isocomma
...
```

---

## 2) Software requirements

The launcher is designed for an HPC module environment and loads:

- BCFtools
- VCFtools
- R
- Biopython

It also requires paths to be set inside `Relate_launcher_V4` for:

- SHAPEIT4
- RELATE
- TWISST

---

## 3) Parameters

### 3.1 Mandatory parameters

- `-v <vcf>`  
  Input VCF (`.vcf.gz`)

- `-c <chromosome>`  
  Contig, scaffold, or chromosome name

- `-s <start_position>`  
  Start coordinate in the VCF

- `-e <end_position>`  
  End coordinate in the VCF

- `--snps <list_snps>`  
  Comma separated SNP positions for which trees are generated

- `-f <pop_file>`  
  Population file (see **1.2**)

- `-t <list_taxa>`  
  Comma separated list of taxa to include in the analysis

- `-o <path_output>`  
  Output folder path

- `-n <name>`  
  Prefix for output files (output folder is `<path_output>/<name>`)

### 3.2 Optional parameters

- `-r <root_outgroup>`  
  Comma separated taxon label(s) used to define the ancestral allele.  
  If not provided, the script uses all samples to define the ancestral state.

- `--species <species_list>`  
  Comma separated list of species labels used as fallback to define the ancestral allele when outgroup data are missing at a SNP.  
  Required when `-r` is provided.

- `--introg <two_species_names>`  
  Comma separated names of two taxa used for the TWISST permutation introgression test.

- `--ps <peak_start>` and `--pe <peak_end>`  
  Peak boundaries used by `twisst_permutations.R`.

---

## 4) Outputs

All outputs are written to:

```text
<output_path>/<name>/
```

Typical retained outputs:

- `*.pdf` (tree plots and TWISST plots)
- `*.newick` (Newick trees)
- `*.pos` (Relate position files)
- `Phylogeny.weights.tsv.gz` (TWISST weights)
- `Other_files/` (archived intermediate outputs)

---

## 5) Pipeline steps (internal logic)

### Step 1: Create ingroup phenotype file

The script creates:

- `<name>_phenotype_file.txt`

by selecting all samples matching the taxa listed in `-t` from the population file.

The script requires at least **20** ingroup samples. If fewer than 20 samples are found, the pipeline exits.

### Step 2: Subset the VCF to the region and ingroup samples

The script subsets the VCF to:

- region: `CHR:START-END`
- samples: from `<name>_phenotype_file.txt`

Filters applied:

- SNPs only
- biallelic only (`-m2 -M2 -v snps`)
- missingness filter: `F_MISSING < 0.2`

Outputs:

- `<name>.vcf.gz` and `<name>.vcf.gz.tbi`

### Step 3: Subset outgroup VCF (only if `-r` is provided)

If `-r` is provided, the script:

1. creates `<name>_outgroup_phenotype_file.txt`
2. subsets the same region and keeps the same SNP positions as the ingroup VCF
3. writes `<name>_outgroups.vcf.gz`

If `-r` is not provided, the script uses all samples to define the ancestral state.

### Step 4: Phase the ingroup VCF (SHAPEIT4)

The script phases:

- input: `<name>.vcf.gz`
- output: `<name>_phased.vcf.gz` and index

### Step 5: Convert phased VCF to RELATE inputs

1. Rename chromosome to `1` (RELATE expects numeric chromosomes)
2. Create a uniform recombination map using `createuniformrecmap.r`
3. Convert to RELATE input format using `RelateFileFormats --mode ConvertFromVcf`

Outputs:

- `<name>.haps`
- `<name>.sample`

### Step 6: Infer ancestral state and flip alleles in the `.haps`

#### Step 6.1 Choose VCF used for ancestral state inference

- If `-r` is provided: use `<name>_outgroups.vcf.gz`
- Otherwise: build `<name>_ALL_samples.vcf.gz` from all samples in the region

#### Step 6.2 Estimate allele frequencies (VCFtools)

The script runs `vcftools --freq` on the chosen ancestral state VCF.

In outgroup mode, it builds a dynamic `--indv` list from outgroup samples.

#### Step 6.3 Assign ancestral allele per SNP

From `*.frq`, the script assigns:

- if REF frequency is missing (`-nan`) → `nan`
- if REF frequency > 0.5 → REF is ancestral
- if REF frequency == 0.5 → `nan`
- else ALT is ancestral

Output:

- `ancestral.alleles`

The script then filters `ancestral.alleles` to keep only SNP positions present in the `.haps`.

#### Step 6.4 Fallback ancestral state for SNPs marked `nan` (species based)

For each species label in `--species`, the script:

1. selects up to **3 random individuals** for that species from the population file
2. subsets the VCF to those individuals and to SNPs present in the ingroup `.vcf.gz`
3. computes allele frequencies (`vcftools --freq`)
4. keeps only SNPs that are `nan` in the main ancestral allele file and are present in the `.haps`

Finally, the script runs:

- `Find_Alternative_ancestral.r`

to update `ancestral.alleles` for sites where outgroup inference was missing.

#### Step 6.5 Flip alleles and recode genotypes in `.haps`

The script creates:

- `<name>_ancestral_state.haps`

For each SNP, if the ancestral allele matches the ALT allele in the `.haps`, it:

- swaps REF and ALT
- flips genotype codes (0 ↔ 1) across all haplotypes

### Step 7: Create RELATE poplabels file

The script creates:

- `<name>_relate.poplabels`

with header:

```text
sample	population	group	sex
```

and assigns sex = `0` for all samples.

### Step 8: Run RELATE

RELATE is executed in `<output_path>/<name>/` using:

- `--mode All`
- mutation rate `-m 2.9e-9`
- effective population size `-N 20000000`
- input: `<name>_ancestral_state.haps`, `<name>.sample`, and the map file

Outputs include:

- `<name>.anc`
- `<name>.mut`
- additional RELATE outputs produced by `--mode All`

### Step 9: Extract trees and plot SNP specific trees

#### Step 9.1 Determine extraction interval

From `--snps`, the script computes:

- `min_value` (minimum SNP position)
- `max_value` (maximum SNP position)

#### Step 9.2 Extract Newick trees across the interval

The script runs `RelateExtract --mode AncToNewick` to create:

- `<name>_<CHR>_<min>-<max>.newick`
- `<name>_<CHR>_<min>-<max>.pos`

#### Step 9.3 Rename sample labels

The script runs `rename_tree.py` to rename tips using an `ID_species` mapping derived from:

- `<name>_phenotype_file.txt`

#### Step 9.4 Plot trees for each SNP

For each SNP position in `--snps`, the script:

1. runs RELATE `TreeView.sh` for the SNP
2. extracts the SNP tree with `tree_extract.py`
3. replots with `plot_tree.R`

Outputs are prefix based and include per SNP PDFs.

### Step 10: Run TWISST

#### Step 10.1 Build position intervals

The script creates `positions.txt` (start and end per interval) using:

- `position_4_twisst.py`

#### Step 10.2 Build TWISST group options

The launcher auto generates `-g` options from unique group names in `group.txt`.

#### Step 10.3 Run TWISST

The script runs TWISST and writes:

- `Phylogeny.weights.tsv.gz`

### Step 11: Plot TWISST results and optional introgression test

If `--introg` is provided, the script:

- parses the two taxa
- runs `twisst_permutations.R` using:
  - `START`, `END`
  - `--ps`, `--pe` peak boundaries

If `--introg` is not provided, it runs `twisst.R`.

### Step 12: Cleanup intermediate files

The script removes intermediate files (VCFs, maps, temporary ancestral state files), then moves remaining non final outputs into:

- `Other_files/`

---

## 6) Notes and constraints

- Output directory `<output_path>/<name>` is deleted and recreated if it already exists.
- If `-r` is provided, `--species` must also be provided.
- Missing SNPs in the ingroup `.haps` are excluded from ancestral state inference.
- The pipeline requires at least **20** ingroup samples to proceed.
