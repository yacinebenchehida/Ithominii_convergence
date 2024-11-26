# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to read and clean data
read_and_clean_data <- function(file_path) {
  data <- read.table(file_path, header = TRUE)
  return(data)
}

# Capture external arguments for the file paths
args <- commandArgs(trailingOnly = TRUE)

# Check if two arguments are provided (for the two input files)
if(length(args) != 2) {
  stop("Please provide two file paths as arguments: file1 (NCD data), file2 (NCDmid data)")
}

# Read the two input data files using the provided arguments
file1 <- args[1]  # First file: NCD, NCDopt, NCDsub
file2 <- args[2]  # Second file: NCDmid

data1 <- read_and_clean_data(file1)  # First dataset with NCD, NCDopt, NCDsub
data2 <- read_and_clean_data(file2)  # Second dataset with NCDmid

# Check for columns in data1 and data2 to ensure they match the expected structure
print("Columns in data1:")
print(names(data1))
print("Columns in data2:")
print(names(data2))

# Prepare the first dataset for plotting (focus on midpos and NCD values)
data1_long <- data1 %>%
  select(midpos, HKA, NCD, NCDopt, NCDsub) %>%
  gather(key = "Statistic", value = "value", -midpos)

# Prepare the second dataset for plotting (focus on snpPos and NCDmid)
data2_long <- data2 %>%
  select(snpPos, NCDmid) %>%
  rename(midpos = snpPos)  # Rename snpPos to match the naming convention in data1

# Create a new row for NCDmid in data1_long
data2_long <- data2_long %>%
  mutate(Statistic = "NCDmid") %>%
  select(midpos, Statistic, value = NCDmid)  # rename NCDmid to value for consistency


# Plot the data
p1 <- ggplot(data1_long, aes(x = midpos, y = value)) +
  geom_point(aes(color = Statistic), size = 2) +  # Add points to the plot
  theme_bw() +
  labs(x = "Mid Position", y = "Statistic", title = "NCD Plot across Mid Positions") +
  facet_wrap(~Statistic, ncol = 1) +  # Separate by NCD types
  theme(legend.position = "none")

p2 <- ggplot(data2_long, aes(x = midpos, y = value)) +
  geom_point(aes(color = Statistic), size = 2) +  # Add points to the plot
  theme_bw() +
  labs(x = "Mid Position", y = "Statistic", title = "NCD Plot across Mid Positions") +
  theme(legend.position = "none")

# Print the plot
pdf("HKA_NCD_plots.pdf",10,8)
plot(p1)
plot(p2)
dev.off()
