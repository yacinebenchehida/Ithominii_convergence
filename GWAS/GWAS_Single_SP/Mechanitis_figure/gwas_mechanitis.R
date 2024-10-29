#!/usr/bin/env Rscript

##############################
# Install required libraries #
##############################
library(ggfortify)    
library(qqman)        
library(tidyverse)    
library(ggtext)       
library(normentR)     
library(grid)         
library(gridExtra)    
library(dplyr)        
library("sjlabelled") 
library(ggplot2)      
library(cowplot)      
library(quotidieR)    
library("ggthemes")   
library("gggenes")    
library(patchwork)    
library(png)          
library(viridis)      
library(wesanderson)  
library(magrittr)     

#################################
# Define a few useful functions #
#################################

# Prepare_data function processes GWAS data for analysis
Prepare_data <- function(gwas) {
  # Rename columns to standard GWAS format
  colnames(gwas) = c("CHR","SNP","BP","P")
  
  # Convert chromosome column to a factor
  gwas$CHR = factor(gwas$CHR)
  
  # Subset significant and non-significant data based on p-value
  sig_data <- gwas %>% subset(P < 1)
  notsig_data <- gwas %>% subset(P >= 1) %>% group_by(CHR)
  
  # Combine both subsets
  gwas_data <- bind_rows(sig_data, notsig_data)

  # Calculate cumulative base pair positions for each chromosome
  data_cum <- gwas_data %>%
    group_by(CHR) %>%
    summarise(max_bp = max(BP)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    select(CHR, bp_add)

  # Merge cumulative positions back into the data
  gwas_data <- gwas_data %>%
    inner_join(data_cum, by = "CHR") %>%
    mutate(bp_cum = BP + bp_add)

  # Set axis positions for each chromosome based on cumulative base pairs
  axis_set <- gwas_data %>%
    group_by(CHR) %>%
    summarize(center = mean(bp_cum))

  # Determine y-axis limit for plots based on minimum p-value
  ylim <- gwas_data %>%
    filter(P == min(P)) %>%
    mutate(ylim = abs(floor(log10(P))) + 2) %>%
    pull(ylim)

  print("data ready")
  
  return(list(gwas_data = gwas_data, axis_set = axis_set))
}

# find_peak function identifies the SNP with the lowest p-value (the GWAS peak)
find_peak <- function(gwas_data){
  thres_pval = 8.108249  # Set threshold p-value
  # Extract peak chromosome and position
  peak = as.character(gwas_data[gwas_data$P == min(unlist(gwas_data$P)), 1])
  peak_position = gwas_data[gwas_data$P == min(unlist(gwas_data$P)), 3]
  peak = peak[1]  
  peak_position = peak_position[1]
  results <- list(peak = peak, peak_position = peak_position)
  return(results)
}

# assess_fixed_mutations checks for genotype consistency across phenotypes
assess_fixed_mutations <- function(my_df){
  unique_phenos <- unique(my_df$pheno)
  
  # Three distinct phenotypes check
  if (length(unique_phenos) == 3){
    genotype1 <- unique(my_df[my_df$pheno == 2, 4])
    genotype2 <- unique(my_df[my_df$pheno == 3, 4])
    genotype3 <- unique(my_df[my_df$pheno == 2.5, 4])
    
    # Ensure unique genotypes for each phenotype
    if (length(genotype1) == 1 & length(genotype2) == 1 & length(genotype3) == 1) {
      if (genotype1 != genotype2 & genotype2 != genotype3 & genotype1 != genotype3) {
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else{
      return(FALSE)
    }
  # Two distinct phenotypes check
  }else if (length(unique_phenos) == 2){
    genotype1 <- unique(my_df[my_df$pheno == 2, 4])
    genotype2 <- unique(my_df[my_df$pheno == 3, 4])
    
    # Check for no overlap in genotypes
    if (length(genotype1) < 3 & length(genotype2) < 3){
      if (!any(genotype1 %in% genotype2)){ 
        print("there are fixed mutations youpi")
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# Fixed_mutations function identifies SNPs within the peak interval
Fixed_mutations <- function(gwas_data, bonf_threshold, peak, peak_position, VCF, phenotype) {
  # Extract significant SNPs and their positions
  signi_SNPs = gwas_data[-log10(gwas_data$P) > 8.108249,]
  write.table(signi_SNPs$BP, file = "Mechanitis_messenoides_tmp.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)
  
  # Define start and end of peak interval
  start_peak <- head(signi_SNPs, n = 1)
  end_peak <- tail(signi_SNPs, n = 1)
  peak_interval <- list(start_peak = start_peak[, 3], end_peak = end_peak[, 3])

  # Run external script for fixed mutations
  command <- paste("./get_fixed_mutations.sh", VCF, peak, peak_interval$start_peak, peak_interval$end_peak, i, sep = " ")
  system(command)
  
  # Read genotype-phenotype file for the phenotype
  geno_pheno <- read.table(paste(phenotype, "_genotype_phenotype_input.txt", sep = ""))
  colnames(geno_pheno) <- c("SNP", "SAMPLE", "pheno", "geno")
  list = list()
  a = 1

  # Loop through each SNP to assess for fixed mutations
  for (i in unique(geno_pheno$SNP)){
    current_snp <- geno_pheno[geno_pheno$SNP == i,]
    current_snp <- current_snp[current_snp$pheno != -9,]
    current_snp$geno <- gsub("/|\\|", "", current_snp$geno)
    current_snp$geno <- gsub("10", "01", current_snp$geno)
    pct_missing <- 100 * (dim(current_snp[current_snp$geno == "..",])[1] / dim(current_snp[current_snp$geno != "..",])[1])
    current_snp <- current_snp[current_snp$geno != "..",]
  
    # Append SNP if it meets criteria for fixed mutations
    if(assess_fixed_mutations(current_snp) == TRUE){
      list[[a]] = cbind(i)
      a = a + 1
    }
  }
  
  list_of_fixed_mutations <- do.call(rbind, list)
  return(list_of_fixed_mutations)
  print(paste("the fixed mutations are: ", list_of_fixed_mutations, sep = ""))
}

# Function to calculate association measures for each SNP
Spear_Cram <- function(phenotype){ 
  geno_pheno <- read.table(paste(phenotype, "_genotype_phenotype_input.txt", sep = ""))
  colnames(geno_pheno) <- c("SNP", "SAMPLE", "pheno", "geno")
  list = list()
  a = 1

  # Loop through each SNP
  for (i in unique(geno_pheno$SNP)){
    current_snp <- geno_pheno[geno_pheno$SNP == i,]
    current_snp <- current_snp[current_snp$pheno != -9,]
    current_snp$geno <- gsub("/|\\|", "", current_snp$geno)
    current_snp$geno <- gsub("10", "01", current_snp$geno)
    current_snp <- current_snp[current_snp$geno != "..",]
    current_snp$GenotypeNumeric <- as.numeric(factor(current_snp$geno, levels = c("00", "01", "11")))
    
    # Calculate Spearman and Cramer's V
    Spearman <- cor(current_snp$GenotypeNumeric, as.numeric(current_snp$pheno), method = "spearman")
    contingency_table <- table(current_snp$geno, current_snp$pheno)
    chi2_result <- chisq.test(contingency_table)
    cramers_v <- sqrt(chi2_result$statistic / (sum(contingency_table) * (min(nrow(contingency_table), ncol(contingency_table)) - 1)))
    list[[a]] = c(i, Spearman^2, cramers_v)
    a = a + 1
  }  

  association <- do.call(rbind, list)
  colnames(association) <- c("Position", "Spearman", "CramerV")
  return(association)
}

# Adds a scale bar to the GWAS plot
add_scale_2_plot <- function(p, data, size){
  scale_bar_length <- size  # Define scale bar size in kb
  difference <- (max(data$Position) - min(data$Position)) * 0.15
  x_position <- (min(data$Position) + difference) / 1000
  print(paste("scale bar position is: ", x_position, sep = ""))
  y_position <- max(-log10(data$P)) * 0.97
  
  # Add scale bar annotation
  scale_bar <- annotate("segment", x = x_position - scale_bar_length, xend = x_position,
                        y = y_position, yend = y_position,
                        color = "black", size = 1)
  text_annotation <- annotate("text", x = x_position - scale_bar_length * 0.5, y = y_position - 2, label = paste0(scale_bar_length, " kb"), size = 5)
  scale_plot <- p + scale_bar + text_annotation
  return(scale_plot)
}

###############
###############
# Main script #
###############
###############
# Load raw GWAS input data
gwas = read.table("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_mechanitis/Data/MechData.txt")

# Initialize lists to store GWAS data for each Mechanitis messenoides phenotype
list_gwas <- list()
list_spearman <- list()
a = 1

# Loop over phenotypes
for (i in c("ForewingBase", "ForewingTip", "Hindwing")){
  # Prepare raw GWAS data
  current_gwas = gwas[gwas$V5 == i, -c(5)]
  ready_data <- Prepare_data(current_gwas)
  gwas_data <- ready_data$gwas_data
  # Add phenotype information to each gwas
  gwas_data$phenotype = i
  # Identify peak position
  peak_info <- find_peak(gwas_data)
  peak = peak_info$peak
  peak_position = peak_info$peak_position
  # Load VCF information for genotype information
  VCF = "/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/Mechanitis_messenoides/Mechanitis_chromosome_5.vcf.gz"
  
  # Find fixed mutations and calculate association measures
  # Test if there are fixed mutations
  fixed_mutations = Fixed_mutations(gwas_data, bonf_threshold, peak, peak_position, VCF, i)
  # Test the strength of the association between phenotype and genotype using the spearman coefficient (and Cramer's V)
  association <- Spear_Cram(i)
  association <- as.data.frame(association)
  association$Phenotype = i
  
  # Filter significant SNPs for plotting
  # Extract SNPs from the GWAS data that are above the threshold of significance (going to be colorised by spearman coefficient in the plot)
  spearman_snps <- as.numeric(unlist(association$Position))
  tmporarement <- gwas_data[gwas_data$BP %in% spearman_snps, 4]
  association$P = tmporarement
  
  list_gwas[[a]] = gwas_data
  list_spearman[[a]] = association  
  a = a + 1  
}

# Combine all data into data frames for plotting
association <- as.data.frame(do.call(rbind, list_spearman))
raw_dt <- as.data.frame(do.call(rbind, list_gwas))
raw_dt <- raw_dt[, -c(5, 6)]
colnames(raw_dt) <- c("Scaff", "whatever", "Position", "P", "Phenotype")


# Generate GWAS plot with color scaling based on Spearman correlation
p <- ggplot(raw_dt, aes(x = Position / 1000, y = -log10(P))) +
  geom_point(color = "gray") +
  scale_x_continuous() +
  labs(x = NULL, y = NULL) +
  facet_grid(Phenotype ~ ., scales = "free_y") +
  geom_hline(yintercept = 8, linetype = "dashed", color = "orange", size = 0.8) + 
  theme_bw() +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 12),
        text = element_text(size = 14),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),  
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(colour = "black", size = 1.5),
        strip.text = element_blank()) +
  geom_point(data = association, aes(x = Position / 1000, y = -log10(P), color = Spearman)) + 
  scale_color_gradientn(colors = c("blue", "gold", "red", "black"), limits = c(0.3, 1)) + 
  geom_hline(yintercept = 8.108249, linetype = "dashed", color = "orange", size = 0.8) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(panel.background = element_rect(colour = "black", size = 1.5))

# Save the plot to PDF with scale bar
pdf("Mechanitis_gwass.pdf", 6, 8)  
plot(add_scale_2_plot(p, raw_dt, 20)) 
dev.off()
