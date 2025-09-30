library(admixtools)
library(tidyverse)

# Args:
# 1: subspecies file (list of populations for this test)
# 2: results directory
# 3: STAT ("f3" or "f4")

args <- commandArgs(trailingOnly = TRUE)
subspecies_file <- args[1]
results_dir <- args[2]
stat_type <- args[3]

# Read populations
pops <- read.table(subspecies_file, stringsAsFactors = FALSE)
pop1 <- pops[1,1]
print(pop1)
pop2 <- pops[2,1]
print(pop2)
pop3 <- pops[3,1]
print(pop3)
pop4 <- ifelse(nrow(pops) >= 4, pops[4,1], NA)

# Prepare PLINK input prefix
data_plink <- paste0(results_dir, "/Inputs")

# Determine which statistic to compute
if(stat_type == "f4") {
    cat("Running f4 statistics...\n")
    out <- f4(data_plink, pop1, pop2, pop3, pop4)
    outfile <- paste0(results_dir, "/", pop1,"_",pop2,"_",pop3,"_",pop4,"_f4.txt")
} else if(stat_type == "f3") {
    cat("Running f3 statistics...\n")
    out <- qp3pop(data_plink, pop1, pop2, pop3)
    outfile <- paste0(results_dir, "/", pop1,"_",pop2,"_",pop3,"_f3.txt")
} else {
    stop("STAT must be either 'f3' or 'f4'")
}

# Write results
write_tsv(out, outfile)
cat("Done. Results saved to", outfile, "\n")
