# Load required library
library(ggplot2)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
input_file <- args[1]            # First argument: path to input file
start_position <- as.numeric(args[2]) # Second argument: start of red square region
end_position <- as.numeric(args[3])   # Third argument: end of red square region
output_name <- paste(args[4],"/nucleotide_diversity.pdf",sep="")

# Load the data from the input file
data <- read.table(input_file, skip=1)

# Create the plot
plot <- ggplot(data, aes(x = V1, y = V2)) + theme_minimal() +
  geom_point() + # Scatter plot of gene position vs. nucleotide diversity
#  geom_line() +
  annotate(
    "rect",
    xmin = start_position, xmax = end_position,
    ymin = -Inf, ymax = Inf,
    fill = "red", alpha = 0.4  # Transparency set directly here
  ) + ylim(0, NA) +
  labs(title = "Nucleotide Diversity Across Gene Positions",
       x = "Gene Position", y = "Nucleotide Diversity") 

# Print the plot
pdf(output_name)
print(plot)
dev.off()
