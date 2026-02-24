#!/usr/bin/env Rscript
###############################################################################
# 02_add_PCA_to_QTL.R
#
# Adds PCA scores (from 01_patternize_PCA.R) as phenotype rows to an R/qtl csvr file.
#
# Input:
#   - QTL file (csv, no header)
#   - PCA scores CSV with SpecimenID + PC columns
#   - Optional ID mapping file to translate photo IDs -> QTL IDs
#
# Output:
#   - 02_QTL_with_PCA.csv  (used by 03_run_QTL_scan.R)
###############################################################################

suppressPackageStartupMessages({
  library(tools)
})

# --- USER SETTINGS ------------------------------------------------------------
qtl_file   <- "dat.4.qtl.lumped.csv"
pca_file   <- "01_PCA_scores.csv"
out_file   <- "02_QTL_with_PCA.csv"

n_pheno_rows <- 12  # number of phenotype rows at top of QTL file

# Optional mapping file (set NULL to skip)
id_map_file <- "match_photo_IDs_QTLs_HW.csv"
map_from_col <- "scanned_ID_only"
map_to_col   <- "vlookup"
drop_values  <- c("#N/A", "n/a", "no_genotype_data")

missing_fill <- "-"

# --- Read inputs --------------------------------------------------------------
if (!file.exists(qtl_file)) stop("QTL file not found: ", qtl_file)
if (!file.exists(pca_file)) stop("PCA file not found: ", pca_file)

qtl_dat <- read.csv(qtl_file, header = FALSE, stringsAsFactors = FALSE)
pca_dat <- read.csv(pca_file, header = TRUE, stringsAsFactors = FALSE)

if (!"SpecimenID" %in% names(pca_dat)) {
  if ("X" %in% names(pca_dat)) {
    names(pca_dat)[names(pca_dat) == "X"] <- "SpecimenID"
  } else {
    stop("PCA file must contain 'SpecimenID' (or legacy 'X').")
  }
}

# --- Optional mapping ---------------------------------------------------------
if (!is.null(id_map_file)) {
  if (!file.exists(id_map_file)) stop("ID mapping file not found: ", id_map_file)

  id_map <- read.csv(id_map_file, header = TRUE, stringsAsFactors = FALSE)
  if (!all(c(map_from_col, map_to_col) %in% names(id_map))) {
    stop("Mapping file must contain columns: ", paste(c(map_from_col, map_to_col), collapse = ", "))
  }

  id_map[[map_to_col]][id_map[[map_to_col]] %in% c("#N/A", "n/a")] <- "no_genotype_data"
  lookup <- setNames(id_map[[map_to_col]], id_map[[map_from_col]])

  mapped <- unname(lookup[pca_dat$SpecimenID])
  has_map <- !is.na(mapped)
  pca_dat$SpecimenID[has_map] <- mapped[has_map]

  pca_dat <- pca_dat[!(pca_dat$SpecimenID %in% drop_values) & !is.na(pca_dat$SpecimenID), , drop = FALSE]
}

# --- Split QTL blocks ---------------------------------------------------------
if (nrow(qtl_dat) < (n_pheno_rows + 1)) stop("n_pheno_rows too large for QTL file.")

qtl_pheno <- qtl_dat[seq_len(n_pheno_rows), , drop = FALSE]
qtl_geno  <- qtl_dat[(n_pheno_rows + 1):nrow(qtl_dat), , drop = FALSE]
qtl_ids   <- as.character(qtl_dat[2, ])

# --- Align PCA to QTL ID order ------------------------------------------------
pc_cols <- setdiff(names(pca_dat), "SpecimenID")
if (length(pc_cols) == 0) stop("No PC columns found in PCA file.")

aligned <- pca_dat[match(qtl_ids, pca_dat$SpecimenID), c("SpecimenID", pc_cols), drop = FALSE]
aligned$SpecimenID <- NULL
aligned[is.na(aligned)] <- missing_fill
aligned <- as.data.frame(lapply(aligned, as.character), stringsAsFactors = FALSE)

# --- Convert to phenotype-row format -----------------------------------------
pc_block <- t(as.matrix(aligned))
pc_block <- as.data.frame(pc_block, stringsAsFactors = FALSE)
pc_block$phenotype <- rownames(pc_block)
pc_block <- pc_block[, c("phenotype", setdiff(names(pc_block), "phenotype")), drop = FALSE]
colnames(pc_block) <- colnames(qtl_pheno)

qtl_pheno_new <- rbind(qtl_pheno, pc_block)

# Legacy formatting: blank columns 2:3 if present (safe, optional)
if (ncol(qtl_pheno_new) >= 3) qtl_pheno_new[, 2:3] <- ""

qtl_new <- rbind(qtl_pheno_new, qtl_geno)

write.table(
  qtl_new,
  file      = out_file,
  sep       = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

message("Done. Wrote: ", out_file)
