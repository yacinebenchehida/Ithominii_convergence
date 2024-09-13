##################
# Load libraries #
##################
library(ComplexHeatmap)
library(colorRamp2)
library("wesanderson")
library("magrittr")
library("dplyr")
library(getopt)
library(viridis)

####################################
# Define the options specification #
####################################
spec <- matrix(c(
  'species',         's', 1, 'character',  'Species name',
  'gene',            'g', 1, 'character',  'Gene name',
  'gwas',            'w', 1, 'character',  'GWAS file',
  'vcf_focal',       'f', 1, 'character',  'VCF file for focal species',
  'vcf_multi',       'm', 1, 'character',  'VCF file for multi species',
  'scaffold',        'r', 1, 'character',  'Scaffold name',
  'start_pos',       'x', 1, 'integer',    'Start position',
  'end_pos',         'y', 1, 'integer',    'End position',
  'phen_focal',      'p', 1, 'character',  'Phenotype file for focal species',
  'phen_multi',      'q', 1, 'character',  'Phenotype file for multi species',
  'threshold',       't', 1, 'numeric',    'Threshold value',
  'species_order',   'o', 1, 'character',  'order file with the samples in the heatmap',
  'color_phenotypes','c', 1, 'character',  'color file: two columns file giving the phenotype and associated color' 
), byrow = TRUE, ncol = 5)


#####################
# Parse the options #
#####################
opt <- getopt(spec)


##############################################
# Check if all required options are provided #
##############################################
if (is.null(opt$species) || is.null(opt$gene) || is.null(opt$gwas) || 
    is.null(opt$vcf_focal) || is.null(opt$vcf_multi) || is.null(opt$scaffold) || 
    is.null(opt$start_pos) || is.null(opt$end_pos) || 
    is.null(opt$phen_focal) || is.null(opt$phen_multi) || 
    is.null(opt$threshold) || is.null(opt$species_order) || is.null(opt$color_phenotypes)){
  cat("Usage: Rscript ./heatmap.R -s <species> -g <gene> -w <gwas> -f <vcf_focal> -m <vcf_multi> -r <region|scaffold|chromosome> -x <start_pos> -y <end_pos> -p <phen_focal> -q <phen_multi> -t <threshold> -o <samples_order_file> -c <color_phenotypes>\n")
  quit(status=1)
}

#################
# Print options #
#################
cat("Species name: ", opt$species, "\n")
cat("Gene name: ", opt$gene, "\n")
cat("GWAS file: ", opt$gwas, "\n")
cat("VCF file for focal species: ", opt$vcf_focal, "\n")
cat("VCF file for multi species: ", opt$vcf_multi, "\n")
cat("Scaffold name: ", opt$scaffold, "\n")
cat("Start position: ", opt$start_pos, "\n")
cat("End position: ", opt$end_pos, "\n")
cat("Phenotype file for focal species: ", opt$phen_focal, "\n")
cat("Phenotype file for multi species: ", opt$phen_multi, "\n")
cat("Threshold: ", opt$threshold, "\n")
cat("Samples order file: ", opt$species_order, "\n")
cat("Color file: ", opt$color_phenotypes, "\n")

####################################
# Set a couple of useful functions #
#################################### 
# Function to get phenotype - genotype for focal and multisp data
phenotype_genotype <- function(species, gene, gwas, vcf_focal, vcf_multi, scaffold, start_pos, end_pos, phen_focal, phen_multi){
	command <- paste("./get_phenotype_genotype_file.sh",species, gene, gwas, vcf_focal, vcf_multi, scaffold, start_pos, end_pos, phen_focal, phen_multi,sep=" ")
	system(command)
}

# Function to prepare focal data
prepare_focal_data <- function(input,threshold){
	data <- read.table(input)
	data = data[data$V3 != -9, ]
	colnames(data) = c("SNP", "Sample_name", "Subspecies", "Phenotype","Genotype","pvalue")
	data$Genotype = gsub(pattern = "/", x = data$Genotype, replacement = "", perl = TRUE)
	data$Genotype = gsub(pattern = "\\|", x = data$Genotype, replacement = "", perl = TRUE)
	data[data$Genotype == "10", 4] = "01"
	data$pvalue = -log10(data$pvalue)
	data <- data[order(data$SNP, data$Subspecies), ]
	data <- data[data$pvalue > threshold, ]
	focal_data <- data

	return(focal_data)	
}


# Function to prepare multisp data
prepare_multisp_data <- function(input,threshold,sp){
	data_multisp <- read.table(input)
	colnames(data_multisp) <- c("SNP","Sample_name", "Subspecies", "Phenotype","Genotype","pvalue")
	
	# Remove VCF separators of genotypes (i.e. / and |)
	data_multisp$Genotype = gsub(pattern = "/", x = data_multisp$Genotype, replacement = "", perl = TRUE)
	data_multisp$Genotype = gsub(pattern = "\\|", x = data_multisp$Genotype, replacement = "", perl = TRUE)

	#Take care of the situation where there are multiallelic states SNPs (Any non biallelic states is called other)
	data_multisp[!(data_multisp$Genotype %in% c("..", "00", "01", "11")),5] = "Other"

	# Extract SNPs from the multisp VCF that are present in the GWAS analysis 
	data_multisp$pvalue = -log10(data_multisp$pvalue)
	data_multisp  <- data_multisp[data_multisp$pvalue > threshold, ]

	# Remove focal species samples from the multi sp data
	species <- gsub(pattern = "\\w+_(\\w+)", x = sp, replacement = "\\1", perl = TRUE)
	data_multisp <- data_multisp[!grepl(species, data_multisp[[3]]), ]

	# Order sample per group so they appear sorted by group in the heatmap
	a = 1
	list = list()
	for (i in unique(data_multisp$SNP)){
	  tmpor <- data_multisp[data_multisp$SNP==i,]
	  tmpor <- tmpor[order(tmpor$Phenotype),]
	  list[[a]] = tmpor
	  a = a + 1
	}
	
	data_multisp <- as.data.frame(do.call(rbind,list))

	return(data_multisp)
}


# Function to merge multisp and focal data
merge_focal_multisp <- function(data, data_multisp, ordre) {
    list = list()
    counter = 1
	
    # Merge data
    for (i in unique(data$SNP)) {
        df1 <- data[data$SNP == i, ]
        df2 <- data_multisp[data_multisp$SNP == i, ]
        tmp <- rbind(df1, df2)
        list[[counter]] <- tmp
        counter = counter + 1
    }

    merged_gwas_multisp <- as.data.frame(do.call(rbind, list))

    # Check if any Subspecies in merged_gwas_multisp is missing from ordre
    ordre_species <- read.table(ordre)
    ordre_species <- as.vector(ordre_species$V1)
	
    missing_samples <- setdiff(as.character(unlist(unique(merged_gwas_multisp$Subspecies))), as.character(unlist(ordre_species)))

    if (length(missing_samples) > 0) {
        warning("Warning: The following samples are missing from the order list:\n", 
                paste(missing_samples, collapse = ", "))
    }

    # Arrange the data as per the order in SNP and Subspecies
    merged_gwas_multisp <- merged_gwas_multisp %>%
        arrange(match(SNP, unique(data$SNP)), match(Subspecies, ordre_species))

    return(merged_gwas_multisp)
}


# Function to create the matrix input of complexheatmap
create_matrix_input <- function(merged_gwas_multisp, ordre){		
	ordre_species <- read.table(ordre, stringsAsFactors = FALSE)
	ordre_species <- as.vector(ordre_species$V1)
	list_snps = list()
	counter1 = 1
	counter2 = 1
	for (i in unique(merged_gwas_multisp$SNP)){
		list_group = list()
		for(j in ordre_species){
		tmp <- merged_gwas_multisp[merged_gwas_multisp$SNP == i & merged_gwas_multisp$Subspecies==j, ]
		list_group[[counter1]] = tmp[,c(5)]
		counter1 = counter1 + 1
	 }
	 list_snps[[counter2]] <- do.call(c,list_group)
	 counter2 = counter2 + 1
	}
	
	Input_matrix = as.matrix(do.call(cbind, list_snps))

	return(Input_matrix)
}

# Function to create the heatmap annotation
create_heatmap <- function(merged_gwas_multisp,colors_file,ordre,Input){
	ordre_species <- read.table(ordre)
	ordre_species <- as.vector(ordre_species$V1)
	# 1) Subspecies annotation
	a = table(merged_gwas_multisp$Subspecies) / length(unique(merged_gwas_multisp$SNP))
	a  <- a[match(unique(unlist(merged_gwas_multisp$Subspecies)), names(a))]
	a = a[a != 0]
	list = list()
	list2 = list()
	counter = 1
	color_sp = as.data.frame(cbind(as.character(unique(merged_gwas_multisp$Subspecies)),rainbow(length(unique(merged_gwas_multisp$Subspecies))))) # These colors are not displayed in the final heatmap
	colnames(color_sp) = c("Species","Color")
	couleur = color_sp[color_sp$V1 %in% unique(merged_gwas_multisp$Subspecies),2]
	for (i in 1:length(a)) {
		list[[counter]] = rep(names(a[i]), a[i])
		tmp = color_sp[color_sp$Species==names(a[i]),2]
		names(tmp) = color_sp[color_sp$Species==names(a[i]),1]
		list2[[counter]] = tmp
		counter = counter + 1
	}
	current_group_order <- do.call(c,list)	
	factor_levels <- factor(current_group_order, levels = ordre_species)
	reordered_group <- factor_levels[order(factor_levels)]
	
	print("Subspecies annotation ready")
	
	# 2) pvalues annotation
	lespvalues <- unique(merged_gwas_multisp$pvalue)
	color_function <- colorRamp2(seq(0,max(unique(merged_gwas_multisp$pvalue))+2,3), rev(magma(length(seq(0,max(unique(merged_gwas_multisp$pvalue))+2,3)))))
	
	list3 <- list()
	counter = 1 

	for (i in unique(merged_gwas_multisp$SNP)){
		lespvalues <- unique(merged_gwas_multisp[merged_gwas_multisp$SNP==i,6])
		list3[[counter]] = lespvalues
		counter = counter + 1
	}
	pval <- do.call(rbind,list3)
	
	list4=list()
	counter = 1
	for (i in 1:length(pval)){
		list4[[counter]]= pval[i]
		counter = counter + 1
	}
	
	print("pvalues annotation ready")

	# 3) phenotype annotation
	a = table(merged_gwas_multisp$Phenotype) / length(unique(merged_gwas_multisp$SNP))
	a = a[a != 0]
	list5 = list()
	list6 = list()
	counter = 1
	colors_corresp <- read.table(colors_file)
	sub <- merged_gwas_multisp[merged_gwas_multisp$SNP==unique(merged_gwas_multisp$SNP)[1],3]
	
	for (i in 1:length(sub)) {
		list5[[counter]] = sub[i]
		info_group <-  merged_gwas_multisp[i,4] 
		tmp = colors_corresp[colors_corresp$V1==info_group,2]
		names(tmp) = c(sub[i])
		list6[[counter]] = tmp
		counter = counter + 1
	}
	
	print("Phenotype annotation ready")
		
	# 4) Set genotype colors
	genotypes <- c("..","00", "01", "11", "Other")
	pal <- c("gray75",wes_palette("Zissou1", 3, type = "continuous"),"pink")
	genotype_color_map <- setNames(pal, genotypes)
	
	get_genotype_colors <- function(genotypes, genotype_color_map) {
	return(genotype_color_map[genotypes])
	}
	
	genotypes_data <- unique(merged_gwas_multisp$Genotype)
	genotype_colors <- get_genotype_colors(genotypes_data, genotype_color_map)
	pal = c("gray75",pal[seq(unique(merged_gwas_multisp$Genotype))-1])
	print("Genotype colors set")

	# 5) Add SNP name at the bottom
	
	colnames(Input) <- unique(merged_gwas_multisp$SNP)
	rownames(Input) <- merged_gwas_multisp[merged_gwas_multisp$SNP==unique(merged_gwas_multisp$SNP)[1],3]
	print(Input)

	# 6) Put it all together
	## Set row annotation (phenotypes/species)
	ha_row = rowAnnotation(Group=do.call(c,list5),  border = TRUE, show_legend = FALSE,  col = list(Group = do.call(c,list6)))

	## Set column annotation (SNPs p-values)
	ha_column = HeatmapAnnotation(pvalue = as.numeric(do.call(c,list4)), col = list(pvalue = color_function),border = TRUE,show_annotation_name = FALSE)
	
	## Create heatmap
	brut <- Heatmap(
	Input,
	show_row_dend = FALSE,
	show_column_dend = FALSE,
	name = "zygosity",
	show_row_names = FALSE,
	show_column_names = TRUE,
	use_raster = FALSE,
	row_title_gp = gpar(fontsize = 8),
	row_split = reordered_group,
	col = genotype_colors,
	border = TRUE, left_annotation = ha_row,
	top_annotation = ha_column, column_names_rot = 45, row_title_rot = 0, 
	column_names_side = "top") 
	
	lgd_boxplot = Legend(labels = colors_corresp$V1, title = "Phenotype",legend_gp = gpar(fill = colors_corresp$V2))
	
	# 7) Final heatmap
	pdf(paste("Heatmap_",opt$species, "_",opt$gene,"_zoomed_out_phylo.pdf",sep=""),11,11)
	draw(brut, heatmap_legend_side="right", annotation_legend_side="right",column_gap = unit(2, "mm"),  # Space between column blocks
     row_gap = unit(2, "mm"), heatmap_legend_list = list(lgd_boxplot))
	dev.off()
     
    return(brut)
}

########
# Main #
########
# Create row genotype - phenotype file
phenotype_genotype(opt$species, opt$gene, opt$gwas, opt$vcf_focal, opt$vcf_multi, opt$scaffold, opt$start_pos, opt$end_pos, opt$phen_focal, opt$phen_multi)

# prepare data for focal species 
my_focal_species_input <- paste(opt$species,"_focal_genotype_phenotype_input.txt",sep="")
focal <- prepare_focal_data(my_focal_species_input,opt$threshold)

# Prepare data for multispecies
my_multi_species_input <- paste(opt$species,"_multisp_genotype_phenotype_input.txt",sep="")
multi <- prepare_multisp_data(my_multi_species_input,opt$threshold,opt$species)

# Combine focal and multispecies
Combined <- merge_focal_multisp(focal, multi, opt$species_order)

# Make complexheatmap matrix input
matrix_input <- create_matrix_input(Combined, opt$species_order)

# Generate annotations, heatmap and plot
my_heatmap <- create_heatmap(Combined,opt$color_phenotypes,opt$species_order,matrix_input)

# Usage example:
# Rscript ./heatmap.R -s Melinaea_mothone -g Cortex -w /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/Cortex/GWAS/Melinaea_mothone.txt -f /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/Melinaea_mothone/GWAS.Melmotiso.base.max0.7N.minGQ10.minQ10.GWASInd.mac2.varbi.CHR4.vcf.gz -m /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/5_Filtering/Results/multisp/Melinaea_mothone/*vcf.gz -r SUPER_4 -x 1385004 -y 1398720 -p /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/introgression_heatmaps/Data/Cortex/Melinaea_mothone/Phenotypes.txt -q /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/introgression_heatmaps/Data/Cortex/Melinaea_mothone/Multisp.txt -t 7.5 -o ../Data/Cortex/Melinaea_mothone/order_file.txt -c ../Data/Cortex/Melinaea_mothone/Color.txt
