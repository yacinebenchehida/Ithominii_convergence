# Load required libraries
library(ggtree)
library(ggplot2)
library(ape)  # For reading the Newick file
library(optparse)  # For handling command-line arguments
library(magrittr)

# Define command-line options
option_list <- list(
  make_option(c("-t", "--tree"), type = "character", default = NULL,
              help = "Path to Newick tree file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "tree_plot.pdf",
              help = "Output file for the tree plot [default: %default]", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if tree file is provided
if (is.null(opt$tree)) {
  print_help(opt_parser)
  stop("Error: Please provide a Newick tree file using the -t option.", call. = FALSE)
}

# Function to extract species name after the last underscore "_"
extract_species <- function(label) {
  parts <- strsplit(label, "_")[[1]]
  species <- tail(parts, 1)  # Get the last part after "_"
  return(species)
}

# Load Newick tree from the provided file
tree <- read.tree(opt$tree)

# Extract species from each tip label
species <- sapply(tree$tip.label, extract_species)

# Create a data frame with tip labels and corresponding species
tree_data <- data.frame(label = tree$tip.label, species = species, stringsAsFactors = FALSE)
tree_data$label <- gsub("Sample_", "", tree_data$label)            # Remove "Sample_" prefix
tree_data$label <- sub("_[^_]+$", "", tree_data$label)  
tree$tip.label <- gsub("Sample_", "", tree$tip.label)     
tree$tip.label <- sub("_[^_]+$", "", tree$tip.label)  

# Plot the tree using ggtree and map the species information directly
p <- ggtree(tree) %<+% tree_data + 
  geom_tiplab(aes(label = label, color = species), align = TRUE, size = 3) +  # Align labels and color by species
  scale_color_manual(values = rainbow(length(unique(tree_data$species)))) +
 theme(legend.position = c(0.15, 0.85)) +
  coord_cartesian(clip = "off")  +
  xlim(0, max(tree$edge.length) * 2.5)

# Save the plot to a file (PDF by default)
ggsave(opt$output, plot = p)

cat("Tree plot saved to", opt$output, "\n")
