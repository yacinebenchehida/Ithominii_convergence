# Load the 'tidyr' library for data manipulation functions
library(tidyr)

# 'args' is used to capture command-line arguments. 
# The first argument is a list of species separated by commas, which is split into a vector.
# The second argument is the path to the data files.
args <- commandArgs(trailingOnly = TRUE)
species_list <- strsplit(args[1], ",")[[1]]  # Split species into a list
path_2_data <- args[2]  # Store the path to data

# Function to combine data from species-specific files
Combine_data <- function(species_list, path_2_data) {
  list = list()  # Initialize an empty list to store data from different species
  sp = 1  # Index for the list
  
  # Loop through each species in the species_list
  for (species in species_list) {
    # Find the species-specific file(s) with ".frq" extension and store the full file path
    species_frequencies_file <- list.files(path = path_2_data, pattern = paste0("^", species, ".*\\.frq$"), full.names = TRUE)
    # Read the file into a dataframe
    frequencies <- read.table(species_frequencies_file)
    # Separate the 'V5' column (e.g., Allele1:frequency1) into two columns: Allele1 and frequency1
    frequencies <- separate(frequencies, V5, into = c("Allele1", "frequency1"), sep = ":")
    # Separate the 'V6' column (e.g., Allele2:frequency2) into two columns: Allele2 and frequency2
    frequencies <- separate(frequencies, V6, into = c("Allele2", "frequency2"), sep = ":")
    # Add a new column to track the species name
    frequencies$Species = species
    # Remove the 'V3' and 'V4' columns from the dataframe as they are not needed
    frequencies <- frequencies[,-c(3,4)]
    # Filter out rows where 'frequency1' has the value '-nan'
    frequencies <- frequencies[frequencies$frequency1 != "-nan",]
    # Store the cleaned data for this species in the list
    list[[sp]] = frequencies
    sp = sp + 1  # Increment the index
  }
  
  # Combine all dataframes in the list into one large dataframe using rbind (row binding)
  complete_df <- as.data.frame(do.call(rbind, list))
  
  # Return the combined dataframe
  return(complete_df)
}

# Function to determine the majority allele based on frequency for each SNP
Get_majority <- function(complete_df) {
  list = list()  # Initialize an empty list to store majority allele information
  SNP = 1  # Index for the list
  
  # Loop through each unique SNP position (V2 column in complete_df)
  for (i in unique(complete_df$V2)) {  
    # Extract all rows corresponding to the current SNP position
    current_SNP <- complete_df[complete_df$V2 == i,]
    # Calculate the sum of frequencies for Allele1 and Allele2 (ignore NA values)
    weight_allele1 <- sum(as.numeric(current_SNP$frequency1), na.rm = TRUE)
    weight_allele2 <- sum(as.numeric(current_SNP$frequency2), na.rm = TRUE)
    # If both allele frequencies are equal, randomly pick one of the alleles
    if (weight_allele1 == weight_allele2) {
      result = c(i, sample(c(unique(current_SNP$Allele1), unique(current_SNP$Allele2)), 1))
    } else {
      # Otherwise, choose the allele with the higher total frequency
      result <- c(i, ifelse(weight_allele1 > weight_allele2, unique(current_SNP$Allele1), unique(current_SNP$Allele2)))
    }
    
    # Store the SNP position and the majority allele in the list
    list[[SNP]] <- result
    SNP = SNP + 1  # Increment the index
  }
  # Combine all results into a dataframe
  SNPs_info <- as.data.frame(do.call(rbind, list))
  # Return the dataframe containing SNP positions and majority alleles
  return(SNPs_info)
}

Set_alternative_ancestral_nan <- function(path_2_data, Alternative_ancestral_state) {
    current_ancestral <- read.table(paste(path_2_data, "/ancestral.alleles", sep=""))
    nan_SNPs <- c()  # Initialize a vector to store SNPs with 'nan'
    
    for (i in 1:nrow(current_ancestral)) {
        current_position <- current_ancestral[i,]
        
        # Check if the ancestral allele is 'nan'
        if (current_position[,2] == "nan") {
            current_SNP <- current_position[,1]
            
            # Retrieve the alternative allele for the current SNP
            alternative_state <- Alternative_ancestral_state[Alternative_ancestral_state$V1 == current_SNP, 2]
            
            # If alternative state is not available, add the SNP to nan_SNPs
            if (length(alternative_state) == 0 || is.na(alternative_state)) {
                current_ancestral[i, 2] <- "nan"
                nan_SNPs <- c(nan_SNPs, current_SNP)  # Append SNP with nan
            } else {
                current_ancestral[i, 2] <- alternative_state
            }
        }
    }
    
    # Return both the updated ancestral alleles and the list of SNPs with 'nan'
    return(list(updated_ancestral = current_ancestral, nan_SNPs = nan_SNPs))
}

State_absent_from_all_taxa_used_4_ancestral <- function(hap_file,Alternative_ancestral_state){
	hap_file <- read.table(hap_file)
	hap_file <- hap_file[,c(3,4,5)]
	Final_Alternative_ancestral_state <- Alternative_ancestral_state$updated_ancestral
	for (i in Alternative_ancestral_state$nan_SNPs){
		Allele1 <- hap_file[hap_file$V3==i,2]
		Allele2 <- hap_file[hap_file$V3==i,3]
		random_Allele <- sample(c(Allele1,Allele2),1)
		Final_Alternative_ancestral_state[Final_Alternative_ancestral_state$V1==i,2] = random_Allele
	}
return(Final_Alternative_ancestral_state[,2])
}

# Combine data for all species
Merged_dt <- Combine_data(species_list, path_2_data)

# Determine the majority ancestral state
Alternative_ncestral_state <- Get_majority(Merged_dt)

# Set alternative ancestral alleles where 'nan' exists
Alternative_ancestral_state <- Set_alternative_ancestral_nan(path_2_data, Alternative_ncestral_state)

# Handle the case were all individuals used as ancestral state are missing data by picking one of the two SNPs from the hap file randomly
Final_ancestral_state <- State_absent_from_all_taxa_used_4_ancestral(args[3],Alternative_ancestral_state)

# Print the new ancestral.state file 
write.table(x = Final_ancestral_state, file = paste(path_2_data,"/ancestral.alleles",sep=""),quote=FALSE,row.names=FALSE, col.names=FALSE)

