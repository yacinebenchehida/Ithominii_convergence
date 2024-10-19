#!/usr/bin/env Rscript
# Script adapted from: https://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/topology-weighting/ and https://github.com/simonhmartin/twisst/blob/master/example_plot.R

# List of functions used to plot twisst results. Copied and pasted from https://github.com/simonhmartin/twisst/blob/master/plot_twisst.R for conveniency.
simple.loess.predict <- function(x, y, span, new_x=NULL, weights = NULL, max = NULL, min = NULL, family=NULL){
    y.loess <- loess(y ~ x, span = span, weights = weights, family=family)
    if (is.null(new_x)) {y.predict <- predict(y.loess,x)}
    else {y.predict <- predict(y.loess,new_x)}
    if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
    if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
    y.predict
    }

smooth.df <- function(x, df, span, new_x = NULL, col.names=NULL, weights=NULL, min=NULL, max=NULL, family=NULL){
    if (is.null(new_x)) {smoothed <- df}
    else smoothed = df[1:length(new_x),]
    if (is.null(col.names)){col.names=colnames(df)}
    for (col.name in col.names){
        print(paste("smoothing",col.name))
        smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, new_x = new_x, max = max, min = min, weights = weights, family=family)
        }
    smoothed
    }

smooth.weights <- function(window_positions, weights_dataframe, span, new_positions=NULL, window_sites=NULL){
    weights_smooth <- smooth.df(x=window_positions,df=weights_dataframe,
                                span=span, new_x=new_positions, min=0, max=1, weights=window_sites)

    #return rescaled to sum to 1
    weights_smooth / apply(weights_smooth, 1, sum)
    }


stack <- function(mat){
    upper <- t(apply(mat, 1, cumsum))
    lower <- upper - mat
    list(upper=upper,lower=lower)
    }

interleave <- function(x1,x2){
    output <- vector(length= length(x1) + length(x2))
    output[seq(1,length(output),2)] <- x1
    output[seq(2,length(output),2)] <- x2
    output
    }


sum_df_columns <- function(df, columns_list){
    new_df <- df[,0]
    for (x in 1:length(columns_list)){
        if (length(columns_list[[x]]) > 1) new_df[,x] <- apply(df[,columns_list[[x]]], 1, sum, na.rm=T)
        else new_df[,x] <- df[,columns_list[[x]]]
        if (is.null(names(columns_list)[x]) == FALSE) names(new_df)[x] <- names(columns_list)[x]
        }
    new_df
    }



plot.weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,density=NULL,lwd=1,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                                        ylab="Weighting", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE){
    #get x axis
    x = positions
    #if a two-column matrix is given - plot step-like weights with start and end of each window    
    if (dim(as.matrix(x))[2]==2) {
        x = interleave(positions[,1],positions[,2])
        yreps=2
        }
    else {
        if (is.null(x)==FALSE) x = positions
        else x = 1:nrow(weights_dataframe)
        yreps=1
        }
    
    #set x limits
    if(is.null(xlim)) xlim = c(min(x), max(x))
    
    #if not adding to an old plot, make a new plot
    if (add==FALSE) plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty)
    
    if (stacked == TRUE){
        y_stacked <- stack(weights_dataframe)
        for (n in 1:ncol(weights_dataframe)){
            y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
            y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
            polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], density=density[n], border=NA)
            }
        }
    else{
        for (n in 1:ncol(weights_dataframe)){
            y = rep(weights_dataframe[,n],each=yreps)
            polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA,density=density[n])
            lines(x,y, type = "l", col = line_cols[n],lwd=lwd)
            }
        }
    }


# Function to check if taxa1 and taxa2 are grouped together in a tree
taxa_are_siblings <- function(tree, taxa1, taxa2) {
  # Get the number of tips (leaf nodes)
  num_tips <- length(tree$tip.label)
  
  # Find the index of taxa1 and taxa2 in the tree
  taxa1_idx <- which(tree$tip.label == taxa1)
  taxa2_idx <- which(tree$tip.label == taxa2)
  
  # Traverse the edges of the tree
  for (node in 1:nrow(tree$edge)) {
    # Get parent node
    parent_node <- tree$edge[node, 1]
    
    # Get all children of this parent node
    child_nodes <- tree$edge[tree$edge[, 1] == parent_node, 2]
    
    # If both taxa are in the child nodes of the same parent, they are siblings
    if (taxa1_idx %in% child_nodes && taxa2_idx %in% child_nodes) {
      return(TRUE)
    }
  }
  
  # If no sibling relationship is found, return FALSE
  return(FALSE)
}

# Main function to read trees and check for sibling taxa
main <- function(taxa1, taxa2, trees) {
   
  # Check each tree and print if taxa1 and taxa2 are siblings
  for (i in seq_along(trees)) {
    tree <- trees[[i]]
    if (taxa_are_siblings(tree, taxa1, taxa2)) {
      return(i)
    }
  }
}

# Function to define plotted area
region2plot <- function(peak_start,peak_end){
  mid <- (peak_start + peak_end)/2
  
  return(list(start_plot= mid - 250000, end_plot = mid + 250000))
}

# BLOCK PERMMUTATIONS (no within block shuffling) - block of constant size user defined
# Permutation Function for Block Permutation
permute_blocks <- function(data, region_start, region_end,block_size) {
  
  data <- data %>%
    mutate(block_id = floor((start - min(start)) / block_size))
  
  # Get unique block_ids and shuffle them
  permuted_block_ids <- sample(unique(data$block_id), replace = FALSE)
  
  permuted_data <- data[,c(5,6)]
  df_permuted <- permuted_data %>%
    arrange(match(block_id, permuted_block_ids))
  
  permuted_data <- cbind(data[,c(1,2,3,4)],df_permuted)
  
  # Subset the permuted data to the region of interest
  region_data <- permuted_data %>%
    filter(start >= region_start & end <= region_end)
  
  # Calculate the frequency of windows with introgression_topo ≥ 0.95
  freq_high_weight <- sum(region_data$introgression_topo >= 0.95)
  
  return(freq_high_weight)
}

# Run the block permutation test multiple times (e.g., 10000 permutations) to get p-value and zscore
block_evaluation <- function(n_permutations,block_size,region_start,region_end,df){
  results <- replicate(n_permutations, permute_blocks(df, region_start, region_end,block_size))
  mu <- mean(results)
  std <- sd(results)
  
  observed_freq <- df %>%
    filter(start >= region_start & end <= region_end) %>%
    summarise(freq_high_weight = sum(introgression_topo >= 0.95)) %>%
    pull(freq_high_weight)
  
  if(observed_freq==0){
     zscore <- "No intro Topo > 95% in peak" 
     p_value <-	"No intro Topo > 95% in peak" 
  }else{
  zscore <- (observed_freq - mu)/std
  p_value <- mean(results >= observed_freq)
  }
   return(list(zscore=zscore,pvalue=p_value))
}


# BLOCK PERMMUTATIONS (with block shuffling) - block of constant size user defined
# Permutation Function for Block Permutation
permute_blocks_within_shuffling <- function(data, region_start, region_end,block_size) {
  
  data <- data %>%
    mutate(block_id = floor((start - min(start)) / block_size))
  
  # Get unique block_ids and shuffle them
  permuted_block_ids <- sample(unique(data$block_id), replace = FALSE)
  
  permuted_data <- data[,c(5,6)]
  df_permuted <- permuted_data %>%
    arrange(match(block_id, permuted_block_ids)) %>%
    group_by(block_id) %>% 
    mutate(introgression_topo = sample(introgression_topo)) %>%
    ungroup()
  
  permuted_data <- cbind(data[,c(1,2,3,4)],df_permuted)
  
  # Subset the permuted data to the region of interest
  region_data <- permuted_data %>%
    filter(start >= region_start & end <= region_end)
  
  # Calculate the frequency of windows with introgression_topo ≥ 0.95
  freq_high_weight <- sum(region_data$introgression_topo >= 0.95)
  
  return(freq_high_weight)
}

# Run the block permutation test multiple times (e.g., 10000 permutations) to get p-value and zscore
block_intra_shuffling_evaluation_ <- function(n_permutations,block_size,region_start,region_end,df){
  results <- replicate(n_permutations, permute_blocks_within_shuffling(df, region_start, region_end,block_size))
  mu <- mean(results)
  std <- sd(results)
  
  observed_freq <- df %>%
    filter(start >= region_start & end <= region_end) %>%
    summarise(freq_high_weight = sum(introgression_topo >= 0.95)) %>%
    pull(freq_high_weight)
  
  if(observed_freq==0){
     zscore <- "No intro Topo > 95% in peak"
     p_value <- "No intro Topo > 95% in peak" 
  }else{
  zscore <- (observed_freq - mu)/std
  p_value <- mean(results >= observed_freq)
  }
   return(list(zscore=zscore,pvalue=p_value))
}


###############
# Main script #
###############
# color to use
cols = c(
  "#0075DC", #Blue
  "#2BCE48", #Green
  "#FFA405", #Orpiment
  "#5EF1F2", #Sky
  "#FF5005", #Zinnia
  "#005C31", #Forest
  "#00998F", #Turquoise
  "#FF0010", #Red
  "#9DCC00", #Lime
  "#003380", #Navy
  "#F0A3FF", #Amethyst
  "#740AFF", #Violet
  "#426600", #Quagmire
  "#C20088", #Mallow
  "#94FFB5") #Jade

#Import used libraries
library(ape)
library(magrittr)
library(dplyr)
library(tidyr)

# Set arguments (i.e. start and end of plot, plus level of smoothing)
args <- commandArgs(trailingOnly = TRUE)

#read topologies file
topos = read.tree(file="Phylogeny.topos")

#weights file with a column for each topology
weights_file = 'Phylogeny.weights.tsv.gz'

weights = read.table(weights_file, header = T)

#normalise rows so weights sum to 1
weights = weights / apply(weights, 1, sum)
weights[is.na(weights)] <- 0

tot_weight = sort(apply(weights, 2, mean, na.rm=TRUE),decreasing=T)
Contribution = tot_weight/sum(tot_weight)*100

#Set introgressed topology
introgression_topo <- main(args[4], args[5], topos)
print(introgression_topo)

#Assign color to each topology
fixed_colors <- c("grey50", "grey75")
num_topologies <- ncol(weights)
colors <- rep(NA, num_topologies)
colors[introgression_topo] <- "red"
non_introgression_indices <- setdiff(1:num_topologies, introgression_topo)
colors[non_introgression_indices[1]] <- "grey50"
colors[non_introgression_indices[2]] <- "grey75"

# with group names in 'V1' and color hex codes in 'V2'
group_names <- names(weights)  # Extract the group names from 'weights'

# Create a named vector of colors using the 'corresp' data frame
#color_mapping <- setNames(corresp$V2, corresp$V1)
topo_color_map <- setNames(colors, paste0("topo_", 1:num_topologies))

# Extract ranking of the tree with a weight 
good_trees = as.numeric(as.character(gsub("topo", "", names(tot_weight), ignore.case = FALSE, perl = TRUE,fixed = FALSE, useBytes = FALSE)))

# Map the group names in 'weights' to their corresponding colors
#group_colors <- color_mapping[group_names]

# Plot each topology contribution in a barplot
pdf("Contribution.pdf",10,10)
names(Contribution) = gsub("topo", "topo_", names(Contribution), ignore.case = FALSE, perl = TRUE,fixed = FALSE, useBytes = FALSE)
barplot(Contribution[1:length(topos)], col = topo_color_map[good_trees],las=2)
dev.off()

# Plot each topology ranked by frequency
pdf("Principal_topology.pdf",12,8)
par(mfrow = c(3,5), mar = c(1,1,2,1), xpd=NA)
a = 0
for (n in good_trees[1:length(topos)]){
  a = a + 1
  current_topo_label <- paste0("topo_", n)
  current_color <- topo_color_map[[current_topo_label]]
  plot.phylo(topos[[n]], type = "clad", edge.color = current_color, edge.width = 5, label.offset = .4, cex = 1)
  mtext(side = 3, text = paste0(current_topo_label, " (", round(Contribution[a], digits = 1), " %)"))
}
dev.off()

# Extract coordinate for each window and prepare it in the format expected by twisst
window_data_file = 'positions.txt'
window_data = read.table(window_data_file)
colnames(window_data) = c("start","end")
window_data$mid = (window_data$start + window_data$end) / 2
window_data$sites <- window_data$end - window_data$start + 1

# Plots stacked barplot without smoothing
pdf("twisst_barplot_no_smoothing.pdf",12,8)
par(mfrow =c(1,1), mar = c(3,3,1,1), xpd = FALSE)
# Plot barplot
if(length(args) != 7){
	plot.weights(weights_dataframe=weights, positions=cbind(window_data$start,window_data$end),
             line_cols=NA, lwd= 0, fill_cols= colors,stacked=TRUE,xlim =c(as.numeric(args[1]),as.numeric(args[2])))
}else {
  peak_start <- as.numeric(args[6])
  peak_end <- as.numeric(args[7])
  plot_area <- region2plot(peak_start,peak_end)

  plot.weights(weights_dataframe=weights, positions=cbind(window_data$start,window_data$end),
            line_cols=NA, lwd= 0, fill_cols= colors,stacked=TRUE,xlim =c(as.numeric(plot_area$start_plot),as.numeric(plot_area$end_plot)))
 
 # Add peak positions
  rect(xleft=peak_start, ybottom=0, xright=peak_end, ytop=1,
     col=NULL, border="black",lwd=3)
}
dev.off()


# Plot stacked barplot with smoothing
pdf("twisst_barplot_with_smoothing.pdf",12,8)
par(mfrow =c(1,1), mar = c(3,3,1,1), xpd = FALSE)

if(length(args) != 7){
       weights_smooth = smooth.weights(window_positions=window_data$mid, weights_dataframe = weights,
                                span = as.numeric(args[3]), window_sites=window_data$sites)       
       plot.weights(weights_dataframe=weights_smooth, positions=window_data$mid,
             line_cols=cols, fill_cols=colors,stacked=TRUE,xlim = c(as.numeric(args[1]), as.numeric(args[2])))
}else{
	# Need to define the genomics interval that would be smoothed. To be sure the smoothing is constant between different plots.
	initial <- as.numeric(row.names(window_data[window_data$start >= plot_area$start_plot & window_data$end <= plot_area$end_plot,])[1])
	the_end <- as.numeric(row.names(tail(window_data[window_data$start >= plot_area$start_plot & window_data$end <= plot_area$end_plot,],1)))
	position2besmoothed <- window_data[window_data$start >= plot_area$start_plot & window_data$end <= plot_area$end_plot,]
	tobesmoothed <- weights[c(initial:the_end),]

	weights_smooth = smooth.weights(window_positions=position2besmoothed$mid, weights_dataframe = tobesmoothed,
                                span = as.numeric(args[3]), window_sites=position2besmoothed$sites)

	plot.weights(weights_dataframe=weights_smooth, positions=position2besmoothed$mid,
             line_cols=cols, fill_cols=colors,stacked=TRUE,xlim = c(plot_area$start_plot,plot_area$end_plot))
	# Add peak positions    
	rect(xleft=peak_start, ybottom=0, xright=peak_end, ytop=1,
    col=NULL, border="black",lwd=3)
} 
dev.off()

print("BARPLOTS PLOTTED")

if(length(args) == 7){
	print("STARTING PERMUTATIONS")
	introgression_topo <- weights[,c(introgression_topo)]
	introgression_data <- as.data.frame(cbind(window_data,introgression_topo))

	n_permutations <- 50000
	block_size <- 100000 
	
	print("PERMUTATIONS INTER BLOCKS")
	block_between_permutations_results <- block_evaluation(n_permutations,block_size,peak_start,peak_end,introgression_data)
	print("PERMUTATIONS INTER AND INTRA BLOCKS")
	block_between_within_permutations_results <- block_intra_shuffling_evaluation_(n_permutations,block_size,peak_start,peak_end,introgression_data)
	
	results <- as.data.frame(do.call(cbind,c(block_between_permutations_results,block_between_within_permutations_results)))
	colnames(results) = c("zscore_inter","pvalue_inter","zscore_inter_intra","pvalue_inter_intra")
	
	write.table(file = "significance.txt", x = results,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
}
