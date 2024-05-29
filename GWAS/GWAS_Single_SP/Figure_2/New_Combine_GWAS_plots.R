#!/usr/bin/env Rscript

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
  order_contigs = read.table(paste("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/",gene,"/scaffold/scaffold_order_",sp,".txt",sep=""))
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
Keep_only_big_scaffold <- function(input){
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
find_peak <- function(gwas,gwas_data,gene,sp){
	thres_pval = 0.05/dim(gwas)[1]
	if(gene=="Optix"){
		if(sp=="Hypothyris_anastasia"){
			peak="scaffold_17"
			peak_position = 9383523
			print(peak_position)
			results <- list(peak = peak, peak_position = peak_position)
			print(paste("peak found ", peak, ": ", peak_position, sep=""))
			return(results)
		}
		else if(sp=="Melinaea_menophilus"){
			peak="SUPER_5"
			peak_position = gwas_data[gwas_data$CHR==peak,]
                        peak_position=gwas_data[gwas_data$P==min(unlist(gwas_data$P)),3]
			print(peak_position)
			results <- list(peak = peak, peak_position = peak_position)
                        print(paste("peak found ", peak, ": ", peak_position, sep=""))
                        return(results)
		}else{
		peak = as.character(gwas_data[gwas_data$P==min(unlist(gwas_data$P)),1])
		peak_position = gwas_data[gwas_data$P==min(unlist(gwas_data$P)),3]
		results <- list(peak = peak, peak_position = peak_position)
		print(paste("peak found ", peak, ": ", peak_position, sep=""))
		return(results)
		}
	}else{
	peak = as.character(gwas_data[gwas_data$P==min(unlist(gwas_data$P)),1])
	print(gwas_data[gwas_data$P==min(unlist(gwas_data$P)),1])
	peak_position = gwas_data[gwas_data$P==min(unlist(gwas_data$P)),3]
	results <- list(peak = peak, peak_position = peak_position)
	#peak = names(sort(table(gwas_data[gwas_data$P < thres_pval,1]),decreasing=T)[1])
	print(paste("peak found ", peak, ": ", peak_position, sep=""))
	
	return(results)
	}
}

# Fonctions to get fixed mutations
assess_fixed_mutations <- function(my_df){
  unique_phenos <- unique(my_df$pheno)
  
  if (length(unique_phenos) == 3){
    genotype1 <- unique(my_df[my_df$pheno == 2, 4])
    genotype2 <- unique(my_df[my_df$pheno == 3, 4])
    genotype3 <- unique(my_df[my_df$pheno == 2.5, 4])
    
    if (length(genotype1) == 1 & length(genotype2) == 1 & length(genotype3) == 1) {
      if (genotype1 != genotype2 & genotype2 != genotype3 & genotype1 != genotype3) {
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else{
      return(FALSE)
    }
  }else if (length(unique_phenos) == 2){
    genotype1 <- unique(my_df[my_df$pheno == 2, 4])
    genotype2 <- unique(my_df[my_df$pheno == 3, 4])
    
    if (length(genotype1) < 3 & length(genotype2) < 3){
      if (!any(genotype1 %in% genotype2)){ # No overlap condition
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

Fixed_mutations <- function(gwas_data, bonf_threshold,peak,peak_position,VCF,gene,species) {
	orange_points = gwas_data[gwas_data$CHR==peak & -log10(gwas_data$P) > bonf_threshold & gwas_data$BP > (peak_position - 500000) & gwas_data$BP < (peak_position + 500000),]
	write.table(orange_points$BP, file = paste(species,"_tmp.txt",sep=""), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)
	start_peak <- head(orange_points,n=1)
	end_peak <- tail(orange_points,n=1)
	peak_interval <- list(start_peak = start_peak[,3], end_peak = end_peak[,3])
	command <- paste("./get_fixed_mutations.sh",VCF,peak,peak_interval$start_peak,peak_interval$end_peak,gene,species,sep=" ")
	system(command)
	geno_pheno <- read.table(paste(sp,"_genotype_phenotype_input.txt",sep=""))
	colnames(geno_pheno) <- c("SNP","SAMPLE","pheno","geno")
	list = list()
	a = 1
	
	if(species!="Melinaea_marsaeus"){
                        geno_pheno <- geno_pheno[geno_pheno$pheno!=2.5,]
                }

	for (i in unique(geno_pheno$SNP)){
		current_snp <- geno_pheno[geno_pheno$SNP==i,]
		current_snp <- current_snp[current_snp$pheno!=-9,]
		current_snp$geno <- gsub("/|\\|", "", current_snp$geno)
		current_snp$geno <- gsub("10", "01", current_snp$geno)
		current_snp <- current_snp[current_snp$geno!="..",]
  
		if(fixed_mutations(current_snp)==TRUE){
			list[[a]] = cbind(i)
			a = a + 1
		}
	
	}
	list_of_fixed_mutations <- as.data.frame(do.call(rbind,list))
	print(list_of_fixed_mutations)
	return(list_of_fixed_mutations)
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
               axis.title.y = element_blank(),
               axis.title.x = element_markdown(size=14),
               plot.title = element_text(hjust = 0.5,size=12),
		axis.text.y = element_text(size = 14),
		axis.text.x = element_text(size = 14)) +
         geom_hline(yintercept=bonf_threshold, linetype="dashed", color = "orange", size=0.8) + 
         geom_point(data = gwas_data[gwas_data$CHR==peak & -log10(gwas_data$P) > bonf_threshold & gwas_data$BP > (peak_position - 500000) & gwas_data$BP < (peak_position + 500000),], aes(x= bp_cum/1000000, y = -log10(P)),colour="orange")
	orange_points = gwas_data[gwas_data$CHR==peak & -log10(gwas_data$P) > bonf_threshold & gwas_data$BP > (peak_position - 500000) & gwas_data$BP < (peak_position + 500000),]
	print(paste("start peak"))
	print(head(orange_points,n=1))
	print(paste("end peak"))
	print(tail(orange_points,n=1))
	 
         # + ggtitle(sp)

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
	print(dim(gwas_data)[1])

	# Combine manhattan plot and images
	if(dim(gwas_data)[1] > 1000000){
	ggp_image <- p + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(panel.background=element_rect(colour="black", size = 1.5)) # + # Combine plot & image
	inset_element(p = my_image,
                left = LEFT - 0.07,
                bottom = 0.87,
                right = LEFT - 0.04,
                top = 0.90) +
	inset_element(p = my_image2,
                left = LEFT - 0.07,
                bottom = 0.69,
                right = LEFT - 0.04,
                top = 0.72)
                
	return(ggp_image)
	}else{
	print("non jen e veux pas")
	ggp_image <- p + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(panel.background=element_rect(colour="black", size = 1.5)) # + # Combine plot & image
        inset_element(p = my_image,
               left = LEFT + 0.1,
              bottom = 0.87,
                right = LEFT + 13 ,
                top = 0.90) +
        inset_element(p = my_image2,
                left = LEFT + 0.10 ,
                bottom = 0.69,
                right = LEFT  + 0.13,
                top = 0.72)

        return(ggp_image)
	}
}
#+ theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + 
# Function to get a data frame only for the zoom data
extract_zoom <- function(region,peak,gwas_data,sp,sign){
	colnames(region) = c("Species","Scaffold","Feature","From","To")
	start_plot = min(unlist(region[,c(4,5)]))
	end_plot = max(unlist(region[,c(4,5)]))
	gwas_zoom = gwas_data[gwas_data$CHR==peak & gwas_data$BP > (start_plot-1) & gwas_data$BP < (end_plot+1),]
	if(sign > 1){
		return(gwas_zoom)
	}else{
		print(head(gwas_zoom))
		gwas_zoom =  cbind(gwas_zoom[,c(1:3)],rev(gwas_zoom$P),gwas_zoom[,c(5:6)])
		colnames(gwas_zoom)[4] = c("P")
		print(head(gwas_zoom))
		return(gwas_zoom)
	}
}


# Function to add scale bar to the plot
add_scale_2_plot <- function(p,data,size){
  # Define the size of the scale bar (in kb)
  scale_bar_length <- size

  # Define position of scale
  difference <- (max(data$BP)-min(data$BP))*0.15
  x_position <- (min(data$BP) + difference)/1000
  print(paste("scale bar position is: ", x_position,sep = ""))
  y_position <- max(-log10(data$P))*0.93

  # Add scale bar
  scale_bar <- annotate("segment", x = x_position - scale_bar_length, xend = x_position,
             y = y_position, yend = y_position,
             color = "black", size = 1)

  # Adjust the distance between the scale bar and text annotation
  text_annotation <- annotate("text", x = x_position - scale_bar_length * 0.5, y = y_position - 2, label = paste0(scale_bar_length, " kb"), size = 5)

  scale_plot <- p + scale_bar + text_annotation

  return(scale_plot)

}

# GWAS on zoom
Plotting_gwas_zoom <- function(region,peak_position,peak,gwas_data,bonf_threshold,sp,sign){
	colnames(region) = c("Species","Scaffold","Feature","From","To")
	start_plot = min(unlist(region[,c(4,5)]))
	end_plot = max(unlist(region[,c(4,5)]))
	print(paste(start_plot, peak_position, end_plot, peak, sep = "\t"))
	gwas_zoom = gwas_data[gwas_data$CHR==peak & gwas_data$BP > (start_plot-1) & gwas_data$BP < (end_plot+1),]
	print(head(gwas_zoom))
	species = gsub(pattern = "(\\w)\\w+_(\\w+)",replacement = "\\1. \\2",x = sp)
	formatted_sp <- substitute(atop(bold(italic(underline(x))), atop("", "")), list(x = species))

	if(sign > 0){
	paragon = gwas_zoom[-log10(gwas_zoom$P) > bonf_threshold,]
	p_zoom <- ggplot(gwas_zoom, aes(x = BP/1000, y = -log10(P))) +
        geom_point(color="gray") +
        scale_x_continuous() +
        labs(x = NULL, y = "-log<sub>10</sub>(p)") +
        theme_bw() +
        theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), axis.title.y = element_markdown(size=10), axis.text.y = element_text(size = 14)) +
	theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
	geom_point(data = paragon, aes(x= BP/1000, y = -log10(P)),colour="orange") +
	geom_hline(yintercept=bonf_threshold, linetype="dashed", color = "orange", size=0.8) +
	theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(panel.background=element_rect(colour="black", size = 1.5)) 

	print((min(gwas_zoom$BP) + (max(gwas_zoom$BP) - min(gwas_zoom$BP))*0.90)/1000)
	
	#p_zoom <- add_scale_2_plot(p_zoom,gwas_zoom,20)

	#p_zoom <- p_zoom + draw_label(formatted_sp, color = "black", size = 14, angle = 0,x=(min(gwas_zoom$BP) + (max(gwas_zoom$BP) - min(gwas_zoom$BP))*0.88)/1000,y=max(-log10(gwas_zoom$P))*0.92)
	

	return(p_zoom)
	}else{
	gwas_zoom =  cbind(gwas_zoom[,c(1:2)],rev(gwas_zoom$BP),gwas_zoom[,c(4:6)])
	colnames(gwas_zoom)[3] = c("BP")
	paragon = gwas_zoom[-log10(gwas_zoom$P) > bonf_threshold,]
	p_zoom <- ggplot(gwas_zoom, aes(x = BP/1000, y = -log10(P))) +
        geom_point(color="gray") +
        scale_x_continuous() +
        labs(x = NULL, y = "-log<sub>10</sub>(p)") +
        theme_bw() +
	theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), axis.title.y = element_markdown(size=10), axis.text.y = element_text(size = 14)) +
	theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
	geom_point(data = paragon, aes(x= BP/1000, y = -log10(P)),colour="orange") +
        geom_hline(yintercept=bonf_threshold, linetype="dashed", color = "orange", size=0.8) + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(panel.background=element_rect(colour="black", size = 1.5))
	
	p_zoom <- add_scale_2_plot(p_zoom,gwas_zoom,20)
	#p_zoom <- p_zoom + draw_label(formatted_sp, color = "black", size = 14, angle = 0,x=(min(gwas_zoom$BP) + (max(gwas_zoom$BP) - min(gwas_zoom$BP))*0.88)/1000,y=max(-log10(gwas_zoom$P))*0.92) # Add species namae

        return(p_zoom)
	}
}

#  axis.ticks.y=element_blank(), axis.text.y=element_blank()
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
Plotting_syntheny_text <- function(region_genes,sp,sign,couleurs){
	colnames(region_genes) = c("Species","Scaffold","Feature","From","To")
	dummies <- make_alignment_dummies(region_genes,aes(xmin = From, xmax = To, y = Species, id = Feature),on = region_genes[1,3])
	my_colors = couleurs

	if(sign > 0){
	region_genes$positions = (region_genes$To + region_genes$From)/2
	syntheny_text <-ggplot(region_genes, aes(xmin = From , xmax = To,y = Species, fill = Feature)) +
	geom_gene_arrow(arrowhead_width = grid::unit(0, "mm") , arrowhead_height = grid::unit(0, "mm"))  +
	scale_color_manual(values = c("black")) + scale_fill_manual(values = my_colors) + guides(fill = guide_legend(nrow = 2)) +
	theme(legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank()) + 
	theme(panel.border = element_blank()) +
	theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"),panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="black"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,160,0), plot.margin=unit(c(-13,1,-0.5,1), "cm")) + guides(fill = FALSE)
	#geom_text(aes(label = Feature, x = positions,  y = 0.95),angle = 20, hjust = 0.65, vjust = -1.5, size = 3) + guides(fill = FALSE) 
	
	print("p_syntheny ready")
	return(syntheny_text)
	}else{
	list=list()
	a = 1
	for (i in 1:dim(region_genes)[1]){
		tmp = region_genes[i,]
		tmp$size = tmp[,5]-tmp[,4]
		list[[a]] = tmp
		a = a + 1
		}
	region_genes = do.call(rbind,list)
	print("region genes after change")
	print(region_genes)
	region_genes = cbind(region_genes[,1:3],rev(region_genes[,4]),rev(region_genes[,5]),region_genes[,6])
	colnames(region_genes) = c("Species","Scaffold","Feature","From","To","size")
	region_genes$To = region_genes$From + region_genes$size
	region_genes$positions = (region_genes$To + region_genes$From)/2
	print(region_genes)
	syntheny_text <-ggplot(region_genes, aes(xmin = From , xmax = To,y = Species, fill = Feature)) +
        geom_gene_arrow(arrowhead_width = grid::unit(0, "mm") , arrowhead_height = grid::unit(0, "mm"))  +
        scale_color_manual(values = c("black")) + scale_fill_manual(values = my_colors) + guides(fill = guide_legend(nrow = 2)) +
        theme(legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank()) +
        theme(panel.border = element_blank()) +
        theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"),panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="black"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,160,0), plot.margin=unit(c(-13,1,-0.5,1), "cm")) + guides(fill = FALSE)
        #geom_text(aes(label = Feature, x = positions,  y = 0.95),angle = 20, hjust = 0.65, vjust = -1.5, size = 3) + guides(fill = FALSE)
        
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

# Determine the number of arguments to separate
n <- length(dat) - 1  # All arguments except the last one

# Separate the arguments
especes <- dat[1:n]
gene <- tail(dat, 1)

list = list() #  Create an empty list that will contain all the plots
counter = 1 #  Will count the element in the list

if(gene=="Optix"){
	couleurs=magma(10)[c(1,6,9)]
}else if(gene=="Antennapedia"){
	couleurs=magma(10)[c(3,8)]
}else{
	couleurs=viridis(5)
}


for (sp in especes){
	cat(paste("SPECIES: ",sp, " GENE: ", gene, sep=""),"\n")
	setwd(paste("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/",gene,sep=""))
	gwas = read.table(paste("GWAS/",sp,".txt",sep=""),skip=1)
	bonf_threshold = -log10(0.05/as.numeric(dim(gwas)[1]))

	images = system("readlink -f images/*", intern = TRUE)
	annotation = read.table("annotation_info/annotation.txt",sep="\t")
	print(annotation)
	setwd("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Scripts") 
		
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
	peak_info = find_peak(gwas,gwas_data,gene,sp)
	peak = peak_info$peak
	peak_position = peak_info$peak_position
	
	#######################
        # Get fixed mutations #
        #######################
	VCF=paste("/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/",sp, "/*vcf.gz",sep="")
	Fixed_mutations(gwas_data, bonf_threshold,peak,peak_position,VCF,gene,sp)

	####################
	# Plot per species #
	####################
	# Manhattan plot
	p <- Plotting_gwas(gwas_data,sp,bonf_threshold,peak,peak_position,ready_data$axis_set)

	# Add butterflies images
	ggp_image <- add_images(images,sp, p, gwas_data)
	
	#########################
	# GWAS zoom around peak #
	#########################
	region = annotation[annotation$V1==sp,]
	print(region)	
	sign = as.numeric(region[2,4]) - as.numeric(region[1,4])
	p_zoom <- Plotting_gwas_zoom(region,peak_position,peak,gwas_data,bonf_threshold,sp,sign)
	print("p_zoom ready")

	######################################
        # Add images to the zoom around peak #
        ######################################
	zoomed_data <- extract_zoom(region,peak,gwas_data,sp,sign)
        
	##############
	# Add images #
	##############
	ggp_zoom_image <- add_images(images,sp, p_zoom, zoomed_data)
	print("zoom image ready")

	############
        # Syntheny #
        ############
	p_syntheny <- Plotting_syntheny_text(region,sp,sign,couleurs)
	print("syntheny ready")

	######################################	
	# Plot syntheny and GWAS around peak #
	######################################
	#png(file=paste(sp,"_GWAS.png",sep=""),width=750,height=500)
        #plot(plot_grid(ggp_image))
	#grid.arrange(p_zoom,p_syntheny)
        #dev.off()

	pdf(paste(sp,"_Zoom_GWAS.pdf",sep=""),8,9)
	plot(plot_grid(ggp_zoom_image,p_syntheny, ncol=1,align="hv",axis="l"))
	#grid.arrange(p_zoom,p_syntheny)
	dev.off()
	print("around peak pdf and png plotted")
	
	############################
        # Move to the next species #
        ############################
	list[[counter]] = p
	counter = counter + 1
}

###########################################
# Combine all the GWAS in a single figure #
###########################################
#if(gene=="Cortex"){
#	png(file="Cortex_GWAS.png",width=700,height=450,type="cairo")
#		plot(wrap_plots(list, ncol = 1))
#	dev.off()
#}

#if(gene=="Optix"){
#	png(file="Optix_GWAS.png",width=700,height=1300,type="cairo")
#		plot(wrap_plots(list, ncol = 1))
#	dev.off()
#}

ggsave(
  "Cortex_GWAS.png",
  plot(wrap_plots(list, ncol = 1)),
  width = 7,
  height = 5,
  dpi = 1200
)
