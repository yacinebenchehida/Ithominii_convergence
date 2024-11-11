# Load necessary libraries
library(data.table)  # Load data.table for efficient data manipulation
library(tidyr)       # Load tidyr for data manipulation functions
library(dplyr)
library(stringr)
library(stringi)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
species_list <- strsplit(args[1], ",")[[1]]  # Split species into a list
path_2_data <- args[2]  # Path to data

# Combine data from species-specific files
Combine_data <- function(species_list, path_2_data) {
  list_of_data <- list()  # Initialize a list to store data from species
  
  # Loop through each species in the species_list
  for (species in species_list) {
    # Read species-specific files with ".frq" extension
    species_frequencies_file <- list.files(path = path_2_data, pattern = paste0("^", species, ".*\\.frq$"), full.names = TRUE)
    frequencies <- fread(species_frequencies_file)  # Use fread from data.table for efficiency
    
    # Separate Allele1 and Allele2 frequencies
    frequencies <- separate(frequencies, V5, into = c("Allele1", "frequency1"), sep = ":")
    frequencies <- separate(frequencies, V6, into = c("Allele2", "frequency2"), sep = ":")
    frequencies$Species <- species  # Track species name
    frequencies <- frequencies[!frequencies$frequency1 == "-nan", ]  # Filter out '-nan'
    
    list_of_data[[species]] <- frequencies  # Store cleaned data
  }
  
  # Combine all dataframes into one
  complete_df <- rbindlist(list_of_data)  # Use rbindlist for efficiency
  return(complete_df)
}

# Determine the majority allele based on frequency for each SNP
Get_majority <- function(complete_df) {
  # Convert to data.table for efficiency
  setDT(complete_df)
  
  # Calculate total frequencies for each SNP
  majority_info <- complete_df[, .(
    frequency1 = sum(as.numeric(frequency1), na.rm = TRUE),
    frequency2 = sum(as.numeric(frequency2), na.rm = TRUE),
    Allele1 = unique(Allele1),
    Allele2 = unique(Allele2)
  ), by = V2]

  # Determine majority allele
  majority_info[, Majority_allele := ifelse(frequency1 > frequency2, Allele1, 
                                             ifelse(frequency2 > frequency1, Allele2, 
                                             sample(c(Allele1, Allele2), 1)))]

  return(majority_info[, .(V2, Majority_allele)])  # Return relevant columns
}

Set_alternative_ancestral_nan <- function(path_2_data, Alternative_ancestral_state) {
    current_ancestral <- read.table(paste(path_2_data, "/ancestral.alleles", sep=""))
    nan_SNPs <- c()  # Initialize a vector to store SNPs with 'nan'

    for (i in 1:nrow(current_ancestral)) {
        current_position <- current_ancestral[i, ]
        
        # Check if the ancestral allele is 'nan'
        if (current_position[, 2] == "nan") {  # Ensure this references the correct column
            current_SNP <- current_position[, 1]  # Ensure this references the correct column
            
            # Retrieve the alternative allele for the current SNP
            alternative_state <- Alternative_ancestral_state[Alternative_ancestral_state$V2 == current_SNP, 2]  # Use correct column here
            
            # If alternative state is not available, add the SNP to nan_SNPs
            if (length(alternative_state) == 0 || is.na(alternative_state)) {
                current_ancestral[i, 2] <- "nan"
                nan_SNPs <- c(nan_SNPs, current_SNP)  # Append SNP with nan
            } else {
                current_ancestral[i, 2] <- alternative_state
            }
        }
    }
    
    if (length(nan_SNPs) == 0) {  # Check length instead of is.null
        nan_SNPs <- c("There are no nans anymore youpi")
    }
    
    return(list(updated_ancestral = current_ancestral, nan_SNPs = nan_SNPs))
}
	
# Function to determine the proposition of ancestral allele present in each species
Count_ancestral <- function(species, ancestral_state_df) {
    genotype_file <- read.table(paste0(species, ".geno"))
    
    # Remove unwanted characters from genotype columns
    genotype_file[,-c(1:3)] <- lapply(genotype_file[,-c(1:3)], 
                                      function(x) str_replace_all(as.character(x), "[/|.]", ""))
    
    # Combine columns and count "0" and "1"
    genotype_file <- genotype_file %>%
      rowwise() %>%
      mutate(
        combined = paste(c_across(-c(V1, V2, V3)), collapse = ""),
        count_0 = str_count(combined, "0"),
        count_1 = str_count(combined, "1")
      ) %>%
      ungroup() %>%
      select(V1, V2, count_0, V3, count_1)

    # Initialize a list to store results
    result_list <- list()
    
    # Iterate over rows and determine output based on ancestral state
    for (i in 1:nrow(genotype_file)) {
        if (genotype_file[i, 2] == ancestral_state_df[i, 2]) {
            result_list[[i]] <- c(genotype_file[i, 3], genotype_file[i, 5])
        } else {
            result_list[[i]] <- c(genotype_file[i, 5], genotype_file[i, 3])
        }
    }
    
    # Print the result as a data frame
    counts_per_sp <- as.data.frame(do.call(rbind, result_list))
    colnames(counts_per_sp) <- c(paste("n",species,sep=""),paste("x",species,sep=""))
    return(counts_per_sp)
    
}


final_input_mutebass <- function(species_list,ancestral_state_df) {
	list = list()
	a = 1
	for (species in species_list){
		print(species)
		counts <- Count_ancestral(species,ancestral_state_df)
		list[[a]] <- counts
		a = a + 1
	}
	mutebass_input <- as.data.frame(do.call(cbind,list))
	mutebass_input <- as.data.frame(cbind(ancestral_state_df[,c(1)],mutebass_input))
	
	column_names <- vector()  # Initialize an empty vector

for (i in 1:dim(mutebass_input)[2]) { 
    if (i == 1) { 
        column_names[i] <- "position"
    } else if (i %% 2 == 0) {  # Check if i is even
        column_names[i] <- paste0("n", i / 2)
    } else {  # i is odd
        column_names[i] <- paste0("x", ceiling(i / 2))
    }
}
colnames(mutebass_input) <- column_names

return(as.data.frame(mutebass_input))

}

# Combine data for all species
Merged_dt <- Combine_data(species_list, path_2_data)

# Determine the majority ancestral state
Alternative_ancestral_state <- Get_majority(Merged_dt)

# Set alternative ancestral alleles where 'nan' exists
Alternative_ancestral_result <- Set_alternative_ancestral_nan(path_2_data, Alternative_ancestral_state)

# Print results in data frame
mute_bass_input <- final_input_mutebass(species_list,Alternative_ancestral_result$updated_ancestral)
mute_bass_input[] <- lapply(mute_bass_input, function(column) {
  if (is.list(column)) {
    unlist(column)
  } else {
    column
  }
})


write.table(x = mute_bass_input, file = paste(path_2_data, "/input_mutebass", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
