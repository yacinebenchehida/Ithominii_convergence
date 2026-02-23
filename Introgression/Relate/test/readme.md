RELATE Pipeline (Relate_launcher_V4)
This folder contains the full pipeline used to perform RELATE-based genealogical inference, ancestral state reconstruction, SNP-specific tree extraction, and TWISST topology weighting analyses for introgression testing.
The pipeline is implemented in:
Relate_launcher_V4
It is designed to run on an HPC using sbatch.
1. Overview of the Pipeline
This pipeline performs the following operations:
Extract region-specific VCF for selected taxa
Phase genotypes with SHAPEIT4
Convert VCF to RELATE input format (.haps and .sample)
Infer ancestral alleles:
Using outgroup allele frequencies when available
Using fallback species sampling when missing
Flip derived/ancestral alleles in .haps
Run RELATE genealogical inference
Extract marginal trees for SNPs of interest
Rename tree labels
Run TWISST topology weighting
Optionally perform introgression permutation tests
Clean intermediate files
2. Software Requirements
The following software must be available:
SHAPEIT4
RELATE (v1.2.1)
BCFtools
VCFtools
Biopython
TWISST
R (≥ 4.2)
Python3
The script assumes module loading on HPC:
module load BCFtools
module load R
module load VCFtools
module load Biopython
Paths to:
SHAPEIT4
RELATE
TWISST
must be correctly defined inside the launcher script.
3. Full Command Usage
sbatch ./Relate_launcher_V4 \
-v <VCF> \
-c <CHROM> \
-s <START> \
-e <END> \
-f <POP_FILE> \
-t <TAXA1,TAXA2,TAXA3,TAXA4> \
--snps <POS1,POS2,...> \
-r <OUTGROUP> \
--species <SPECIES_LIST> \
-o <OUTPUT_PATH> \
-n <PREFIX> \
--introg <SP1,SP2> \
--ps <PEAK_START> \
--pe <PEAK_END>
4. Required Inputs
4.1 VCF
Must be bgzipped and indexed
Must contain all taxa listed in -t
Must contain SNPs for full region
4.2 Population File
Tab-delimited file:
sample_name    group    species
First column: sample ID (must match VCF)
Species names must match those provided via --species
5. Pipeline Steps (Internal Logic)
Step 1 — Create Phenotype File
Extract individuals belonging to taxa in -t
Require ≥ 20 individuals
Output:
PREFIX_phenotype_file.txt
Step 2 — Subset VCF
Extract:
Region: CHR:START-END
Individuals from phenotype file
Only biallelic SNPs
Missingness filter: F_MISSING < 0.2
Output:
PREFIX.vcf.gz
Step 3 — Outgroup VCF (Optional)
If -r is provided:
Extract matching SNP positions
Subset to outgroup individuals
Output:
PREFIX_outgroups.vcf.gz
If not provided:
Use all individuals to infer ancestral state
Step 4 — Phasing (SHAPEIT4)
shapeit4 --input PREFIX.vcf.gz --region CHR
Output:
PREFIX_phased.vcf.gz
Step 5 — Convert to RELATE Format
Rename chromosome to numeric
Generate uniform recombination map
Convert VCF to:
PREFIX.haps
PREFIX.sample
Step 6 — Ancestral State Inference
6.1 Primary Method (Outgroup)
Compute allele frequencies using VCFtools
If REF frequency > 0.5 → REF ancestral
If ALT frequency > 0.5 → ALT ancestral
If frequency = 0.5 or missing → assign nan
Output:
ancestral.alleles
Filter to SNPs present in .haps.
6.2 Fallback Method (Species-based)
For SNPs where ancestral state is nan:
For each species in --species:
Select up to 3 random individuals
Estimate allele frequencies
Resolve ancestral state via R script:
Find_Alternative_ancestral.r
6.3 Flip Alleles in .haps
If ancestral allele matches ALT:
Swap REF and ALT
Invert genotype coding (0 ↔ 1)
Output:
PREFIX_ancestral_state.haps
Step 7 — Generate RELATE Poplabels
Creates:
PREFIX_relate.poplabels
Format:
sample   population   group   sex
Sex fixed to 0.
Step 8 — Run RELATE
Parameters:
Mutation rate: 2.9e-9
Effective population size: 20,000,000
Mode: All
Outputs:
PREFIX.anc
PREFIX.mut
PREFIX.coal
...
Step 9 — Extract SNP Trees
Identify min/max SNP positions
Extract marginal trees:
RelateExtract --mode AncToNewick
Output:
PREFIX_CHR_min-max.newick
Extract individual SNP trees
Replot using:
plot_tree.R
Step 10 — TWISST
10.1 Generate position intervals
positions.txt
10.2 Build group option
Auto-generate -g arguments.
10.3 Run TWISST
twisst.py --method complete
Output:
Phylogeny.weights.tsv.gz
Step 11 — Introgression Testing (Optional)
If --introg provided:
Run:
twisst_permutations.R
Using:
Start
End
Peak start
Peak end
If not provided:
Run standard TWISST plotting script.
Step 12 — Cleanup
Removes:
Intermediate VCF files
Temporary ancestral files
Maps
Moves remaining auxiliary files into:
Other_files/
Final outputs retained:
PDF trees
TWISST plots
Newick trees
Significant topology results
6. Scientific Rationale
Each analysis includes:
Two focal taxa (P1, P2)
A third taxon (P3)
A fourth taxon (P4)
Optional outgroup
This mirrors an ABBA-BABA-like framework, where genealogical topologies at GWAS-significant SNPs are evaluated for discordance patterns consistent with introgression.
Marginal trees are inferred for each SNP exceeding significance thresholds.
TWISST is used to quantify topology weighting across the region.
Permutation tests evaluate whether topology enrichment within peak boundaries deviates from expectation.
7. Typical Use Case
Example:
sbatch ./Relate_launcher_V4 \
-v My_VCF.gz \
-c SUPER_4 \
-s 1 \
-e 1843557 \
-f yellow_band_info.txt \
-t mothone_mothone,mothone_messenina,isocomma,simulator \
--snps 1,1390852,1390869,... \
-r outgroup \
--species tarapotensis,satevis,lilis,isocomma,idae,marsaeus,menophilus,mothone,flavo \
-o ../Results \
-n Example_run \
--introg mothone_mothone,isocomma \
--ps 1385004 \
--pe 1398720
8. Important Constraints
VCF must contain all SNPs in focal region.
Missing SNPs in multi-species VCF will cause failure.
At least 20 individuals required in ingroup.
Species names must exactly match population file.
Script overwrites output folder if it exists.
