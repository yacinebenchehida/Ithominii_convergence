#!/usr/bin/env Rscript

# Load necessary libraries
library(ape)
library(phangorn)

# Get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Read the PHYLIP alignment file
alignment <- read.phyDat(args[1], format = "phylip")

# Compute the GTR distance matrix
gtr_model <- dist.ml(alignment, model = "F81")

# Build the NJ tree
nj_tree <- nj(gtr_model)

# Save the tree to a file
path = args[2]
write.tree(nj_tree, file = paste(path,"/nj_tree.newick",sep=""))
