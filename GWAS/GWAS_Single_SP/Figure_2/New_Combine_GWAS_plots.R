##############################
# Install required libraries #
##############################
#install.packages("ggfortify")
library(ggfortify)
#install.packages("qqman")
library(qqman)
library(tidyverse)
#install.packages("ggtext")
library(ggtext)
#devtools::install_github("norment/normentR")
library(normentR)
library(grid)
library(gridExtra)
library(dplyr)
#install.packages("sjlabelled")
library("sjlabelled")
library(ggplot2)
library(cowplot)
#remotes::install_github("HanjoStudy/quotidieR")
library(quotidieR)
#install.packages("ggthemes")
library("ggthemes")
library("gggenes")
library(patchwork)
library(png)
library(viridis)

#################################
# Define a few useful functions #
#################################
# Prepare data for gwas 
Prepare_data <- function(gwas, sp) {
  colnames(gwas) = c("CHR","SNP","BP","P")
  order_contigs = read.table(paste("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/scaffold_order_",sp,".txt",sep=""))
	order_contigs = as.data.frame(order_contigs)
	order_contigs = as.factor(order_contigs$V1)
	gwas$CHR=factor(gwas$CHR, levels=order_contigs)
	sig_data <- gwas %>% subset(P < 1)
	notsig_data <- gwas %>% subset(P >= 1) %>% group_by(CHR) %>% sample_frac(0.01)
	gwas_data <- bind_rows(sig_data, notsig_data)

	data_cum <- gwas_data %>% 
	group_by(CHR) %>% 
	summarise(max_bp = max(BP)) %>% 
	mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
	select(CHR, bp_add)

	gwas_data <- gwas_data %>% 
	inner_join(data_cum, by = "CHR") %>% 
	mutate(bp_cum = BP + bp_add)

	axis_set <- gwas_data %>% 
	group_by(CHR) %>% 
	summarize(center = mean(bp_cum))

	ylim <- gwas_data %>% 
	filter(P == min(P)) %>% 
	mutate(ylim = abs(floor(log10(P))) + 2) %>% 
	pull(ylim)	
	print("data ready")
	
  return(list(gwas_data = gwas_data, axis_set = axis_set))
}

# Function to only keep big scaffolds
Keep_only_big_scaffold <- function(input) {
list = list()
	a = 1
	for (i in unique(input$CHR)){
		Size = length(input[input$CHR==i,1])
		if(Size > 10000){
			list[[a]] = i
			a = a + 1
		}
	}

	Big_scaff = do.call(c,list)
	big_scaff_data = input[input$CHR %in% Big_scaff,]
	print("big scaffold extracted")
	
	return(big_scaff_data)
}

# Function to find gwas peak
find_peak <- function(gwas,gwas_data){
	thres_pval = 0.05/dim(gwas)[1]
	min(unlist(gwas_data$P))
	peak = as.character(gwas_data[gwas_data$P==min(unlist(gwas_data$P)),1])
	peak_position = gwas_data[gwas_data$P==min(unlist(gwas_data$P)),3]
	results <- list(peak = peak, peak_position = peak_position)
	#peak = names(sort(table(gwas_data[gwas_data$P < thres_pval,1]),decreasing=T)[1])
	print(paste("peak found ", peak, ": ", peak_position, sep=""))
	
	return(results)
}

# Function to create the manhattan plot for the whole genome
Plotting_gwas <- function(gwas_data, sp, bonf_threshold, peak, peak_position, axis_set) {
    p <- ggplot(gwas_data, aes(x = bp_cum/1000000, y = -log10(P),color =CHR)) +
         geom_point() +
         scale_x_continuous() +
         scale_color_manual(values = rep(c("gray", "black"), unique(length(axis_set$CHR)))) +
         labs(x = "Genome position (Mb)", y = "-log10(p)") + 
         theme_bw() +
         theme(legend.position = "none",
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank(),
               axis.title.x = element_markdown(size=12),
               axis.title.y = element_markdown(size=12),
               plot.title = element_text(hjust = 0.5,size=12),
               text = element_text(size=14)) +
         geom_hline(yintercept=bonf_threshold, linetype="dashed", color = "orange", size=0.8) + 
         geom_point(data = gwas_data[gwas_data$CHR==peak & -log10(gwas_data$P) > bonf_threshold & gwas_data$BP > (peak_position - 500000) & gwas_data$BP < (peak_position + 500000),], aes(x= bp_cum/1000000, y = -log10(P)),colour="orange") +
         ggtitle(sp)

    return(p)
}

# Function to append images of the butterflies to the manhattan plot
add_images <- function(images, sp, p, gwas_data){
images = images[grepl(sp,images)]
	print(images)
	my_image <- readPNG(images[1], native = TRUE)
	my_image2 <- readPNG(images[2], native = TRUE)
	
	# Set images position in the manhattan plot
	LEFT = (min(gwas_data[order(gwas_data$P)[1:3],6]) - min(gwas_data$bp_cum))/(max(gwas_data$bp_cum) - min(gwas_data$bp_cum))
	RIGHT = (max(gwas_data[order(gwas_data$P)[1:3],6]) - min(gwas_data$bp_cum))/(max(gwas_data$bp_cum) - min(gwas_data$bp_cum))
	print("Images imported")

	# Combine manhattan plot and images
	ggp_image <- p + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +                # Combine plot & image
	inset_element(p = my_image,
                left = LEFT - 0.13,
                bottom = 0.87,
                right = LEFT - 0.10,
                top = 0.90) +
	inset_element(p = my_image2,
                left = LEFT - 0.13,
                bottom = 0.69,
                right = LEFT - 0.10,
                top = 0.72)
                
    return(ggp_image)
}

# GWAS on zoom
Plotting_gwas_zoom <- function(region,peak_position,peak,gwas_data,bonf_threshold,sp){
	colnames(region) = c("Species","Scaffold","Feature","From","To")
	start_plot = min(unlist(region[,c(4,5)]))
	end_plot = max(unlist(region[,c(4,5)]))
	print(paste(start_plot, peak_position, end_plot, peak, sep = "\t"))
	gwas_zoom = gwas_data[gwas_data$CHR==peak & gwas_data$BP > (start_plot-1) & gwas_data$BP < (end_plot+1),]
	
	if(sp!="Hypothyris_anastasia"){
	paragon = gwas_zoom[-log10(gwas_zoom$P) > bonf_threshold,]
	p_zoom <- ggplot(gwas_zoom, aes(x = BP/1000, y = -log10(P))) +
        geom_point(color="gray") +
        scale_x_continuous() +
        labs(x = "Genome position (kb)", y = "-log<sub>10</sub>(p)") +
        theme_bw() +
        theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
	theme(axis.title.x=element_blank()) +
	theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(),plot.margin=unit(c(0.8,1,-1.2,1), "cm")) +
	geom_point(data = paragon, aes(x= BP/1000, y = -log10(P)),colour="orange") +
	geom_hline(yintercept=bonf_threshold, linetype="dashed", color = "orange", size=0.8) 
    
	return(p_zoom)
	}else{
	gwas_zoom =  cbind(gwas_zoom[,c(1:2)],rev(gwas_zoom$BP),gwas_zoom[,c(4:6)])
	colnames(gwas_zoom)[3] = c("BP")
	paragon = gwas_zoom[-log10(gwas_zoom$P) > bonf_threshold,]
	p_zoom <- ggplot(gwas_zoom, aes(x = BP/1000, y = -log10(P))) +
        geom_point(color="gray") +
        scale_x_continuous() +
        labs(x = "Genome position (kb)", y = "-log<sub>10</sub>(p)") +
        theme_bw() +
        theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
        theme(axis.title.x=element_blank()) +
        theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(),plot.margin=unit(c(0.8,1,-1.2,1), "cm")) +
        geom_point(data = paragon, aes(x= BP/1000, y = -log10(P)),colour="orange") +
        geom_hline(yintercept=bonf_threshold, linetype="dashed", color = "orange", size=0.8)

        return(p_zoom)
	}
}

# Plot genes syntheny
Plotting_syntheny <- function(region_genes){
	colnames(region_genes) = c("Species","Scaffold","Feature","From","To")
	dummies <- make_alignment_dummies(region_genes,aes(xmin = From, xmax = To, y = Species, id = Feature),on = region[1,2])
	
	# Set colors
	my_colors = viridis(length(unique(region_genes$Feature)))
	
	# Plot syntheny based on genes only (no features plotted)
	p_syntheny <- ggplot(region_genes, aes(xmin = From , xmax = To,y = Species, fill = Feature)) +
	geom_gene_arrow(arrowhead_width = grid::unit(0, "mm") , arrowhead_height = grid::unit(0, "mm"))  +
	scale_color_manual(values = c("black")) + scale_fill_manual(values = my_colors) + guides(fill = guide_legend(nrow = 2)) +
	theme(legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank()) + 
	theme(panel.border = element_blank()) +
	theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"),panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="black"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,160,0), plot.margin=unit(c(1,1,-0.5,1), "cm")) 
	print("p_syntheny ready")
	
	return(p_syntheny)
}

# Plot genes syntheny with genes names above genes
Plotting_syntheny_text <- function(region_genes,sp){
	colnames(region_genes) = c("Species","Scaffold","Feature","From","To")
	dummies <- make_alignment_dummies(region_genes,aes(xmin = From, xmax = To, y = Species, id = Feature),on = region_genes[1,3])
	my_colors = viridis(length(unique(region_genes$Feature)))
	region_genes$positions = (region_genes$To + region_genes$From)/2 

	if(sp!="Hypothyris_anastasia"){
	syntheny_text <-ggplot(region_genes, aes(xmin = From , xmax = To,y = Species, fill = Feature)) +
	geom_gene_arrow(arrowhead_width = grid::unit(0, "mm") , arrowhead_height = grid::unit(0, "mm"))  +
	scale_color_manual(values = c("black")) + scale_fill_manual(values = my_colors) + guides(fill = guide_legend(nrow = 2)) +
	theme(legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank()) + 
	theme(panel.border = element_blank()) +
	theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"),panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="black"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,160,0), plot.margin=unit(c(-6,1,-0.5,1), "cm")) +
	geom_text(aes(label = Feature, x = positions,  y = 0.925),angle = 20, hjust = 0.65, vjust = -1.5, size = 3) + guides(fill = FALSE) 
	print("p_syntheny ready")
	
	return(syntheny_text)
	}
}

# Plot syntheny based on gene and features
Plotting_syntheny_with_features <- function(region){
	region_genes = region[region$Feature!="stop_codon",]
	region_genes = region_genes[region_genes$Feature!="start_codon",]
	region_genes = region_genes[region_genes$Feature!="gene",]
	dummies <- make_alignment_dummies(region_genes,aes(xmin = Start, xmax = End, y = Species, id = Gene, forward=orientation),on = region_genes[1,2])
	region_genes$alpha = 0.15
	p_syntheny_features <- ggplot(region_genes, aes(xmin = Start , xmax = End,y = Species, forward = orientation, label = Gene)) +
	geom_gene_arrow(aes(color = Gene), fill = "white",arrowhead_width = unit(1.5, "mm"),arrow_body_height= unit(4,"mm")) + scale_color_manual(values = my_colors) +
	geom_subgene_arrow(arrowhead_width = unit(1.5, "mm"),arrow_body_height= unit(4,"mm"),
	data = region_genes,
	aes(xmin = Start, xmax = End, xsubmin = From , xsubmax = To, fill = Feature,forward = orientation), alpha = region_genes$alpha) +
	scale_fill_manual(values = c("grey","white")) +  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1)) +
	theme(legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank()) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "white"),panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="black"),legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,0,160,0), plot.margin=unit(c(1,1,-0.5,1), "cm"))
	print("p_syntheny_features ready")
	
	return(p_syntheny_features)
}

###############
# Import data #
###############
dat = commandArgs(trailingOnly=TRUE)
list = list() #  Create an empty list that will contain all the plots
counter = 1 #  Will count the element in the list

for (sp in dat){
	cat(paste("SPECIES: ",sp,sep=""),"\n")
	gwas = read.table(paste("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/",sp,".txt",sep=""),skip=1)
	bonf_threshold = -log10(0.05/as.numeric(dim(gwas)[1]))
	images = system("readlink -f /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/images/*", intern = TRUE)
	annotation = read.table("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/annotation_info/annotation.txt",sep="\t") 
		
	################          
	# Prepare data #
	################
	ready_data <- Prepare_data(gwas, sp)
	gwas_data <- ready_data$gwas_data
	
	#####################################          
	# Extract scaffold larger than 2 Mb #
	#####################################
	gwas_data <- Keep_only_big_scaffold(gwas_data)
	
	########################################################################
	# Find scaffold with the peak and the associated points above the peak #
	########################################################################
	thres_pval = 0.05/dim(gwas)[1]
	peak_info = find_peak(gwas,gwas_data)
	peak = peak_info$peak
	peak_position = peak_info$peak_position
	print(peak)
	print(peak_position)

	####################
	# Plot per species #
	####################
	# Manhattan plot
	p <- Plotting_gwas(gwas_data,sp,bonf_threshold,peak,peak_position,ready_data$axis_set)
	
	# Add butterflies images
	images = system("readlink -f /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/images/*", intern = TRUE)
	ggp_image <- add_images(images,sp, p, gwas_data)
	
	#########################
	# GWAS zoom around peak #
	#########################
	region = annotation[annotation$V1==sp,]
	print(region)
	p_zoom <- Plotting_gwas_zoom(region,peak_position,peak,gwas_data,bonf_threshold,sp)
	
	############
        # Syntheny #
        ############
	#p_syntheny <- Plotting_syntheny(region)
	p_syntheny <- Plotting_syntheny_text(region,sp)

	######################################	
	# Plot syntheny and GWAS around peak #
	######################################
	png(file=paste(sp,"_Zoom_GWAS.png",sep=""),width=600,height=500)
        grid.arrange(p_zoom,p_syntheny)
        dev.off()

	pdf(paste(sp,"_Zoom_GWAS.pdf",sep=""),9,7)
	grid.arrange(p_zoom,p_syntheny)
	dev.off()
	print("around peak pdf and png plotted")
	
	############################
        # Move to the next species #
        ############################
	list[[counter]] = ggp_image 
	counter = counter + 1
}

###########################################
# Combine all the GWAS in a single figure #
###########################################
png(file="Cortex_GWAS.png",width=600,height=950,type="cairo")
	wrap_plots(list, ncol = 1)
dev.off()
