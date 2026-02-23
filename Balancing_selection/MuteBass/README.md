# MuteBass Pipeline

## Purpose

This pipeline prepares population genomic data for **MuteBaSS** to test
for balancing selection using:

-   HKA
-   NCD
-   NCDopt
-   NCDsub
-   NCDmid

Analyses are performed in sliding windows (1000 bp windows, 500 bp
step).

------------------------------------------------------------------------

## Required Inputs

### 1. VCF file

-   BGZipped VCF (`.vcf.gz`)
-   SNPs only
-   Includes ingroup and outgroup individuals

### 2. Population file

Tab-delimited file:

    sample_ID    species_name

-   One line per individual
-   One species labeled exactly as `outgroup`
-   All other species are treated as ingroup taxa

### 3. Species Tree (Newick format)

Provided via `-t` argument.

-   Taxa must be numbered according to **alphabetical order of species
    names**

-   Example:

        ((((2,3),4),5),((1,6),7))

-   Use `NADA` if no phylogenetic filtering is required

### 4. Target Frequency

Provided via `-f`.

-   Defines the target allele frequency used in NCD statistics

------------------------------------------------------------------------

## Core Workflow

### 1. Subsampling individuals

`select_random.py`\
- Randomly selects **5 individuals per species** - If fewer than 5 are
available, all are retained

This standardizes sampling across taxa.

------------------------------------------------------------------------

### 2. Ancestral State Assignment

-   Outgroup allele frequencies are computed
-   Ancestral allele is defined as:
    -   Major allele if frequency \> 0.5
    -   `nan` if 0.5 or missing
-   When ancestral state is `nan`, the global majority allele across
    species is used as fallback

Output:

    ancestral.alleles

------------------------------------------------------------------------

### 3. Construction of MuteBaSS Input

`Input_MuteBass.r`

For each SNP and each species:

-   `x_i` = count of derived alleles
-   `n_i` = total allele count

Final input format:

    position  x1  n1  x2  n2  x3  n3 ...

Filtering applied:

-   SNPs polymorphic in more than one species are removed
-   Sites with missing data are removed
-   Only SNPs compatible with MuteBaSS assumptions are retained

Output:

    input_mutebass

------------------------------------------------------------------------

### 4. Phylogenetic Compatibility Filtering (Optional)

If a tree is provided:

-   `MuteBaSS.py --check` identifies SNPs incompatible with the species
    tree
-   Incompatible SNPs are removed

Output:

    Final_mutebass_input.txt

If `-t NADA`, no filtering is performed.

------------------------------------------------------------------------

### 5. Sliding Window Analysis

MuteBaSS is run with:

-   Window size: 1000 bp
-   Step size: 500 bp

Statistics computed:

-   **HKA**
-   **NCD**
-   **NCDopt**
-   **NCDsub**
-   **NCDmid**

Target frequency defined with `--tf`.

Primary output:

    final_results.txt

Columns:

    midpos  HKA  NCD  NCDopt  NCDsub  numSites  optF

------------------------------------------------------------------------

### 6. Plotting

`Plotting_results.r`

Generates:

    NCD_plots.pdf

-   One panel per statistic
-   Values plotted across genomic midpoints

------------------------------------------------------------------------

## Final Outputs

Main results:

-   `final_results.txt`
-   `NCD_plots.pdf`

These files contain all statistics used for downstream interpretation.
