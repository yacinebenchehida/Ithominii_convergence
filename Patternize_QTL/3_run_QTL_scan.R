#!/usr/bin/env Rscript
###############################################################################
# 03_run_QTL_scan.R
#
# QTL mapping of phenotype PCs using R/qtl:
#   - read cross from 02_QTL_with_PCA.csv (csvr)
#   - calc genotype probabilities
#   - scanone per PC phenotype column
#   - permutation thresholds
#   - bayesint intervals
#   - write results table + scan plots
###############################################################################

suppressPackageStartupMessages({
  library(qtl)
})

options(bitmapType = "cairo")

# --- USER SETTINGS ------------------------------------------------------------
cross_file <- "02_QTL_with_PCA.csv"

# Genotype coding (adjust for your cross)
genotypes <- c("AA", "AB", "BB")
alleles   <- c("A", "B")

# Genoprob settings
step_cM      <- 1
off_end      <- 0.0
error_prob   <- 0.001
map_function <- "haldane"
stepwidth    <- "fixed"

# Phenotype columns to scan (these must match where PCs ended up in your QTL file)
# If your original file had 12 phenotype rows and you added PC1..PC5, they will be 13:17.
pheno_cols <- 13:17

# Permutations
n_perm <- 1000
alpha_for_threshold <- 0.05

# Bayes interval
bayes_prob <- 0.95

# Marker file for genetic->physical intervals
marker_file <- "linkage_map.csv"
marker_file_has_header <- FALSE

results_dir <- "results"
analysis_label <- "03_QTL_scan"

# --- Helper functions ---------------------------------------------------------
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

safe_chr <- function(x) {
  x <- as.character(x)
  x[x %in% c("21", "chr21")] <- "X"
  x
}

read_marker_table <- function(path, has_header = FALSE) {
  mk <- read.csv(path, header = has_header, stringsAsFactors = FALSE)

  if (all(c("scaff", "pos", "cM", "chr") %in% names(mk))) {
    mk2 <- mk[, c("scaff", "pos", "cM", "chr")]
    mk2$chr <- safe_chr(mk2$chr)
    mk2$cM  <- as.numeric(mk2$cM)
    mk2$pos <- as.numeric(mk2$pos)
    return(mk2)
  }

  # Legacy logic (edit here if your file differs)
  if (!has_header && ncol(mk) >= 7) {
    mk <- mk[8:nrow(mk), c(3, 4, 7)]
    mk$chr <- substr(mk$V3, 6, 7)
    colnames(mk) <- c("scaff", "pos", "cM", "chr")
    mk$chr <- safe_chr(mk$chr)
    mk$cM  <- as.numeric(mk$cM)
    mk$pos <- as.numeric(mk$pos)
    return(mk)
  }

  stop("Marker file format not recognized. Provide scaff,pos,cM,chr or edit read_marker_table().")
}

thin_markers_by_cM <- function(mk) {
  out <- NULL
  for (chr in unique(mk$chr)) {
    chr_dat <- mk[mk$chr == chr, , drop = FALSE]
    for (cm in unique(chr_dat$cM)) {
      cm_dat <- chr_dat[chr_dat$cM == cm, , drop = FALSE]
      if (nrow(cm_dat) == 1) out <- rbind(out, cm_dat)
      if (nrow(cm_dat) > 1) {
        out <- rbind(out, cm_dat[1, ])
        out <- rbind(out, cm_dat[nrow(cm_dat), ])
      }
    }
  }
  out
}

physical_interval_from_bayes <- function(bint, markers_thin) {
  chr_low  <- as.character(bint[1, 1])
  chr_peak <- as.character(bint[2, 1])
  chr_high <- as.character(bint[3, 1])

  cm_low  <- as.numeric(bint[1, 2])
  cm_peak <- as.numeric(bint[2, 2])
  cm_high <- as.numeric(bint[3, 2])

  low_rows  <- markers_thin[markers_thin$chr == chr_low  & round(markers_thin$cM, 3) == round(cm_low, 3), , drop = FALSE]
  peak_rows <- markers_thin[markers_thin$chr == chr_peak & round(markers_thin$cM, 3) == round(cm_peak, 3), , drop = FALSE]
  high_rows <- markers_thin[markers_thin$chr == chr_high & round(markers_thin$cM, 3) == round(cm_high, 3), , drop = FALSE]

  if (nrow(low_rows) == 0 || nrow(peak_rows) == 0 || nrow(high_rows) == 0) {
    return(list(peak_low = NA, peak_high = NA, physical_peak = NA, physical_limits = NA))
  }

  peak_low  <- paste0(peak_rows[1, "scaff"], ":", peak_rows[1, "pos"])
  peak_high <- paste0(peak_rows[nrow(peak_rows), "scaff"], ":", peak_rows[nrow(peak_rows), "pos"])

  peak_pos_vec <- as.numeric(peak_rows$pos)
  median_pos <- peak_pos_vec[ceiling(length(peak_pos_vec) / 2)]
  median_scaff <- as.character(peak_rows$scaff[which(as.numeric(peak_rows$pos) == median_pos)[1]])
  physical_peak <- paste0(median_scaff, ":", median_pos)

  physical_limits <- paste0(
    low_rows[1, "scaff"], ":", low_rows[1, "pos"],
    "-",
    high_rows[nrow(high_rows), "scaff"], ":", high_rows[nrow(high_rows), "pos"]
  )

  list(peak_low = peak_low, peak_high = peak_high, physical_peak = physical_peak, physical_limits = physical_limits)
}

star_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.1)   return("^")
  ""
}

# --- Read cross ---------------------------------------------------------------
if (!file.exists(cross_file)) stop("Cross file not found: ", cross_file)

cross <- read.cross("csvr", file = cross_file, genotypes = genotypes, alleles = alleles)
cross <- jittermap(cross)
cross <- calc.genoprob(
  cross,
  step         = step_cM,
  off.end      = off_end,
  error.prob   = error_prob,
  map.function = map_function,
  stepwidth    = stepwidth
)

if (any(pheno_cols > ncol(cross$pheno))) stop("pheno_cols exceed number of phenotype columns in cross$pheno.")

# --- Markers ------------------------------------------------------------------
if (!file.exists(marker_file)) stop("Marker file not found: ", marker_file)

markers <- read_marker_table(marker_file, has_header = marker_file_has_header)
markers_thin <- thin_markers_by_cM(markers)

# --- Output -------------------------------------------------------------------
results_file <- file.path(results_dir, paste0(analysis_label, "_results.tsv"))
scanplot_prefix <- file.path(results_dir, paste0(analysis_label, "_"))

header <- paste(
  "Cross", "Phenotype",
  "LOD_max", "chrom_max",
  "peak_interval_low", "peak_interval_high",
  "cM_max", "cM_limits",
  "physical_peak", "physical_limits",
  "Intercept", "Term1", "Term2", "Term3",
  "R2",
  sep = "\t"
)
writeLines(header, con = results_file)

# --- QTL scans ----------------------------------------------------------------
for (pheno_col in pheno_cols) {

  pheno_name <- names(cross$pheno)[pheno_col]
  message("Scanning: ", pheno_name)

  scan <- scanone(cross, pheno.col = pheno_col, model = "normal", method = "hk")
  perm <- scanone(cross, pheno.col = pheno_col, n.perm = n_perm, model = "normal", method = "hk", perm.Xsp = FALSE)

  thr_table <- summary(perm, alpha = c(0.10, 0.05, 0.01, 0.001))
  threshold <- as.numeric(thr_table[which(c(0.10, 0.05, 0.01, 0.001) == alpha_for_threshold), 1])

  # Plot scan
  png(filename = paste0(scanplot_prefix, pheno_name, "_scan.png"), width = 1200, height = 400)
  plot(scan, main = paste0("scanone: ", pheno_name))
  abline(h = threshold, col = "red")
  dev.off()

  # Per chromosome
  for (chrom in unique(scan$chr)) {

    temp_chr <- scan[scan$chr == chrom, , drop = FALSE]
    if (nrow(temp_chr) == 0) next

    lod_max_chr <- max(temp_chr$lod, na.rm = TRUE)
    if (!is.finite(lod_max_chr) || lod_max_chr <= threshold) next

    bint <- bayesint(scan, chr = chrom, prob = bayes_prob)
    LOD_max   <- round(bint[2, 3], 2)
    chrom_max <- as.character(bint[2, 1])
    cM_max    <- round(bint[2, 2], 2)
    cM_limits <- paste0(round(bint[1, 2], 2), "-", round(bint[3, 2], 2))

    temp_chr_phys <- temp_chr[!grepl("loc", rownames(temp_chr), fixed = TRUE), , drop = FALSE]
    bint_phys <- bayesint(temp_chr_phys, prob = bayes_prob)
    phys <- physical_interval_from_bayes(bint_phys, markers_thin)

    peak_marker <- find.marker(cross, bint[2, 1], bint[2, 2])

    # Simple LM at peak marker (based on genotype probabilities)
    pheno_vec <- pull.pheno(cross, pheno_col)
    geno_prob <- pull.genoprob(cross, chr = chrom_max, omit.first.prob = FALSE, include.pos.info = TRUE, rotate = TRUE)

    gp_row <- geno_prob[geno_prob$marker == peak_marker, , drop = FALSE]
    if (nrow(gp_row) != 1) next

    gp <- as.data.frame(t(gp_row[, 5:ncol(gp_row), drop = FALSE]))
    if (ncol(gp) >= 2) gp <- gp[, -1, drop = FALSE]  # baseline genotype absorbed into intercept

    model <- lm(pheno_vec ~ ., data = gp)
    sm <- summary(model)
    coefs <- coef(sm)

    coef_strings <- rep(NA_character_, 4) # intercept + up to 3 terms
    names(coef_strings) <- c("Intercept", "Term1", "Term2", "Term3")

    if ("(Intercept)" %in% rownames(coefs)) {
      est <- round(coefs["(Intercept)", "Estimate"], 2)
      se  <- round(coefs["(Intercept)", "Std. Error"], 2)
      p   <- coefs["(Intercept)", "Pr(>|t|)"]
      coef_strings["Intercept"] <- paste0(est, "±", se, star_p(p))
    }

    other_terms <- setdiff(rownames(coefs), "(Intercept)")
    if (length(other_terms) > 0) {
      for (k in seq_len(min(3, length(other_terms)))) {
        term <- other_terms[k]
        est <- round(coefs[term, "Estimate"], 2)
        se  <- round(coefs[term, "Std. Error"], 2)
        p   <- coefs[term, "Pr(>|t|)"]
        coef_strings[paste0("Term", k)] <- paste0(est, "±", se, star_p(p))
      }
    }

    out_line <- paste(
      "cross", pheno_name,
      LOD_max, chrom_max,
      phys$peak_low, phys$peak_high,
      cM_max, cM_limits,
      phys$physical_peak, phys$physical_limits,
      coef_strings["Intercept"], coef_strings["Term1"], coef_strings["Term2"], coef_strings["Term3"],
      round(sm$r.squared, 2),
      sep = "\t"
    )

    write(out_line, file = results_file, append = TRUE)
  }
}

message("Done. Wrote results to: ", results_file)
