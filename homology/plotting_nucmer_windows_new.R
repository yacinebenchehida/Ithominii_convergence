library(dplyr)
library(ggplot2)
library(viridis)


dat = commandArgs(trailingOnly=TRUE)
annot <- read.table(dat[1],sep="\t")
gene <-annot$V3[1]


if(gene=="Hyd. like"){
  Inputs <- c("mapping_nucmer_sliding_windows_Melinaea_marsaeus_Hypothyris_anastasia.txt", "mapping_nucmer_sliding_windows_Melinaea_menophilus_Melinaea_marsaeus.txt",
           "mapping_nucmer_sliding_windows_Mechanitis_messenoides_Melinaea_menophilus.txt")
}else{
  Inputs <- c("mapping_nucmer_sliding_windows_Melinaea_mothone_Hypothyris_anastasia.txt", "mapping_nucmer_sliding_windows_Melinaea_menophilus_Melinaea_mothone.txt",
           "mapping_nucmer_sliding_windows_Mechanitis_messenoides_Melinaea_menophilus.txt")
}



############################################################################
# Function to assess if the Spearman coefficient for each SNPs in the peak #
############################################################################
Spear_square_coeff <- function(sp){
        geno_pheno <- read.table(paste(sp,"_genotype_phenotype_input.txt",sep=""))
        colnames(geno_pheno) <- c("SNP","SAMPLE","pheno","geno")
        list = list()
        a = 1

	for (i in unique(geno_pheno$SNP)){
                current_snp <- geno_pheno[geno_pheno$SNP==i,]
                current_snp <- current_snp[current_snp$pheno!=-9,]
                current_snp$geno <- gsub("/|\\|", "", current_snp$geno)
                current_snp$geno <- gsub("10", "01", current_snp$geno)
                current_snp <- current_snp[current_snp$geno!="..",]
                current_snp$GenotypeNumeric <- as.numeric(factor(current_snp$geno, levels = c("00", "01", "11")))
                Spearman <- cor(current_snp$GenotypeNumeric, as.numeric(current_snp$pheno), method = "spearman")

                list[[a]] = c(i, Spearman^2)
                a = a + 1
        }
	association <- do.call(rbind,list)
        colnames(association) <- c("SNP","Spearman")
        print(association)
        return(association)
}

Spear_square_coeff("Melinaea_menophilus")


#######################################################################
# Function that puts all the data in the right format before plotting #
#######################################################################
prepare_data <- function(aln,sp1, sp2, anno) {

aln = read.table(paste("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Results/mummer/sliding_windows_mapping/",aln,sep=""),header = TRUE,fill = TRUE)
query_species = sp1
target_species = sp2
aln <- aln[aln$Identity > 40,]

###################################
# Load the annotation information #
###################################
anno <- anno[anno$V1 %in% c(target_species,query_species),]

# Function to check the coordinate of the annotation and adjust the value. Also if the reference genome is reverse complemented it reorder the annotation so i$
check_direction_and_realign_annotation <- function(anno){
  if(anno$V4[2] - anno$V4[1] > 0){
    star_position <- min(anno$V4)
    anno$V4 <- anno$V4 - star_position
    anno$V5 <- anno$V5 - star_position
    return(anno)
  }
  else{
    star_position <- min(anno$V4)
    anno$V4 <- anno$V4 - star_position
    anno$V5 <- anno$V5 - star_position
    anno$size <- anno$V5 - anno$V4
    anno$succes_size <- 0
    for (i in 1:(dim(anno)[1]-1)){
      anno$succes_size[i] <- anno$V4[i] - anno$V5[i+1]
    }
    anno$V4[1] = 0
    anno$V5[1] =  anno$size[1]
    for (i in 2:(dim(anno)[1])){
      anno$V4[i] <- anno$V5[i-1] + anno$succes_size[i-1]
      anno$V5[i] <- anno$V4[i] + anno$size[i]
    }
    return(anno[,c(1:5)])
  }
}

#########################################################
# Prepare the properly annotation file for the plotting #
##########################################################
list = list()
a = 1
for (i in c(query_species,target_species)){
  tmp <- anno[anno$V1==i,]
  tmp <- check_direction_and_realign_annotation(tmp)
  #  tmp <- tmp[tmp$V5< max(input$queryEnd),]
  list[[a]] = tmp
  a = a + 1
}

anno = as.data.frame(do.call(rbind,list))

############################
# Figure out genome offset #
############################
for (i in 1:dim(aln)[1]){
  aln[i,9] = aln[i,9] + aln[i,1] - 1 
  aln[i,10] = aln[i,10] + aln[i,1] - 1 
  aln[i,5] = aln[i,5] 
  aln[i,6] = aln[i,6] 
}

##########################################################################################
# Remove interval aligning to very distant non homologous part of the alternative genome #
##########################################################################################                                          
size_difference <- abs(max(anno[anno$V1==target_species,5])-max(anno[anno$V1==query_species,5]))
if(size_difference < 5000){
  size_difference = 5000
}

tolerance <-  seq(from = 7000, to = size_difference +7000, length.out = dim(aln)[1])
list = list()
counter = 1

differential <- max(anno[anno$V1==target_species,5])-max(anno[anno$V1==query_species,5])

for (i in 1:dim(aln)[1]){
  tol <- tolerance[i]
  difference <- aln$subjectEnd[i]-aln$queryEnd[i]
  if(abs(difference) <  tol & differential < 0 & difference < 0){
    print(paste(i,abs(difference), differential, tol))
    list[[counter]] = aln[i,]
    counter = counter + 1
  }
  else if(abs(difference) <  tol & differential > 0 & difference > 0){
    list[[counter]] = aln[i,]
    counter = counter + 1
  }
  else if(abs(difference) < 5000){
    list[[counter]] = aln[i,]
    counter = counter + 1
  }
}

aln <- as.data.frame(do.call(rbind,list))

#########################################################
# Reformat the data in the geom_polygon expected format #
#########################################################
list = list()
counter = 1

for (i in 1:dim(aln)[1]){
  list[[counter]] <- data.frame(
    x = c(aln$subjectStart[i],aln$subjectEnd[i],  aln$queryEnd[i],aln$queryStart[i]),
    y = c(1,1,3,3),
    group = paste("window_",i,sep=""), 
    order = c(1, 2, 4, 3),
    reference = c(target_species, target_species,query_species, query_species))
  counter = counter + 1
}

plotting_data <-  do.call(rbind,list)

#####################################################################################################
# Extract alignment regions in the annotation file. These bits will be colorised in bed in the plot #
#####################################################################################################
to_colorise_query <- anno[anno$V1==query_species,c(3,4,5)]
to_colorise_target <- anno[anno$V1==target_species,c(3,4,5)]

list = list()
counter = 1
for (i in unique(plotting_data$group)){
  tmp1 = plotting_data[plotting_data$group==i & plotting_data$y==1,]
  tmp2 = plotting_data[plotting_data$group==i & plotting_data$y==3,]
  for (j in 1:dim(to_colorise_target)[1]){
    if(any(!is.na(to_colorise_target[j, 2]) & !is.na(to_colorise_target[j, 3]) & tmp1$x >=  to_colorise_target[j, 2] & tmp1$x <= to_colorise_target[j, 3])){
      if(any(!is.na(to_colorise_query[j, 2]) & !is.na(to_colorise_query[j, 3]) & tmp2$x >= to_colorise_query[j, 2] & tmp2$x <= to_colorise_query[j, 3])){
        list[[counter]] = plotting_data[plotting_data$group==i,]
        counter = counter + 1
      }
    }
  }
}


to_colorise <- do.call(rbind,list)

##################
# set color code #
##################
gene <-anno$V3[1]

if(gene=="Hyd. like"){
  color_code <- magma(10)[c(1,6,9,4)]
}else{
  color_code <- c(viridis(5),"orange")
}


return(list(plotting_data, anno, to_colorise, color_code,query_species, target_species))


}


#########################################
# Function to add scale bar to the plot #
#########################################
add_scale_to_plot <- function(p, data, size) {
  # Define the size of the scale bar (in kb)
  scale_bar_length <- size
  
  # Define position of scale
  difference <- (max(data$x) - min(data$x)) * 0.12
  x_position <- min(data$x) + difference
  cat("Scale bar position is:", x_position, "\n")
  y_position <- 8.1
  
  # Add scale bar
  scale_bar <- annotate(
    "segment", 
    x = x_position - scale_bar_length, xend = x_position,
    y = y_position, yend = y_position,
    color = "black", size = 1
  )
  
  # Add text annotation
  text_annotation <- annotate(
    "text", 
    x = x_position - scale_bar_length * 0.5, y = y_position - 0.1, 
    label = paste0(scale_bar_length / 1000, " kb"), size = 5
  )
  
  # Add the annotations to the plot
  scale_plot <- p + scale_bar + text_annotation
  
  return(scale_plot)
}


########################################################
# Define a function to create a layer for each dataset #
########################################################
 create_layer <- function(data_file, y_offset) {
   sp1 <- gsub(x = data_file, pattern = "mapping_nucmer_sliding_windows_(\\w+_\\w+)_(\\w+_\\w+).txt", replacement = "\\1", perl = TRUE)
   sp2 <- gsub(x = data_file, pattern = "mapping_nucmer_sliding_windows_(\\w+_\\w+)_(\\w+_\\w+).txt", replacement = "\\2", perl = TRUE)
   ready_data <- prepare_data(data_file, sp1, sp2, annot)
   plotting_data <- ready_data[[1]]
   anno <- ready_data[[2]]
   to_colorise <- ready_data[[3]]
   color_code <- ready_data[[4]]
   query_species <- ready_data[[5]]
   target_species <- ready_data[[6]]
   
   # Adjust y-axis position for stacking
   plotting_data$y <- plotting_data$y + y_offset
   
   # Create a list of ggplot layers for the current dataset
   layers <- list(
     geom_polygon(data = plotting_data, aes(x = x, y = y, group = group), color = "grey", fill = "black", alpha = 0.4),
     geom_rect(aes(xmin = 0, xmax = max(anno[anno$V1 == target_species, 5]), ymin = 0.95 + y_offset, ymax = 0.75 + y_offset), color = "black", fill = "white"),
     geom_rect(data = anno[anno$V1 == unique(anno$V1)[2], ], aes(xmin = V4, xmax = V5, ymin = 0.948 + y_offset, ymax = 0.752 + y_offset, fill = V3), inherit.aes = FALSE),
     geom_rect(aes(xmin = 0, xmax = max(anno[anno$V1 == query_species, 5]), ymin = 3.05 + y_offset, ymax = 3.25 + y_offset), color = "black", fill = "white"),
     geom_rect(data = anno[anno$V1 == unique(anno$V1)[1], ], aes(xmin = V4, xmax = V5, ymin = 3.052 + y_offset, ymax = 3.248 + y_offset, fill = V3), inherit.aes = FALSE),
     geom_polygon(data = to_colorise, aes(x = x, y = y + y_offset, group = group), color = "brown", fill = "brown")
   )
   
   list(layers = layers, color_code = color_code, y_values = unique(plotting_data$y), species_labels = c(target_species, query_species),plotting_data = plotting_data)
 }
 
 # Apply the function to each dataset, accumulating the y-offset as well
 layers_list <- lapply(1:length(Inputs), function(i) {
   create_layer(Inputs[i], y_offset = (i - 1) * 2.3)
 })
 
 # Extract individual components from the list
 all_layers <- unlist(lapply(layers_list, function(x) x$layers), recursive = FALSE)
 color_code <- layers_list[[1]]$color_code  # Assuming color code is the same across all layers
 y_values <- unlist(lapply(layers_list, function(x) x$y_values))
 species_labels <- unlist(lapply(layers_list, function(x) x$species_labels))
 
 # Combine all layers into a single plot starting from an empty ggplot object
 p <- ggplot() +
   all_layers +
   scale_y_continuous(breaks = y_values[c(1,3,5,6)], labels = species_labels[c(1,3,5,6)]) +
   scale_fill_manual(values = color_code) +
   theme_minimal() +
   labs(x = NULL, y = NULL) +
   theme(axis.text.x = element_blank())


# Get data necessary to add scale 
for_scale <- create_layer(Inputs[1], y_offset = (1 - 1) * 2.3)

# Print the final plot
if(gene=="Hyd. like"){
	pdf("plot_syntheny_windows_nucmer_optix.pdf",12,10)
	plot(add_scale_to_plot(p,for_scale$plotting_data, size = 20000))
	dev.off()
}else{
	pdf("plot_syntheny_windows_nucmer_cortex.pdf",12,10)
	plot(add_scale_to_plot(p, for_scale$plotting_data, size = 20000))
	dev.off()
}
