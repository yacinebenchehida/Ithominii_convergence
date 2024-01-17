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
	
	################          
	# Prepare data #
	################
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
	#####################################          
	# Extract scaffold larger than 2 Mb #
	#####################################
	list2 = list()
	a = 1
	for (i in unique(gwas_data$CHR)){
		Size = length(gwas_data[gwas_data$CHR==i,1])
		if(Size > 10000){
			list2[[a]] = i
			a = a + 1
		}
	}

	Big_scaff = do.call(c,list2)
	gwas_data = gwas_data[gwas_data$CHR %in% Big_scaff,]
	print("big scaffold extracted")
	
	########################################################################
        # Find scaffold with the peak and the associated points above the peak #
        ########################################################################
	thres_pval = 0.05/dim(gwas)[1]
	min(unlist(gwas_data$P))
	peak = gwas_data[gwas_data$P==min(unlist(gwas_data$P)),1]
	peak_position = gwas_data[gwas_data$P==min(unlist(gwas_data$P)),3]
	#peak = names(sort(table(gwas_data[gwas_data$P < thres_pval,1]),decreasing=T)[1])
	print("peak found")

	####################
	# Plot per species #
	####################
	 p <- ggplot(gwas_data, aes(x = bp_cum/1000000, y = -log10(P),color =CHR)) +
	 geom_point() +
	 scale_x_continuous() +
	 scale_color_manual(values = rep(c("gray", "black"), unique(length(axis_set$CHR)))) +
	 labs(x = "Genome position (Mb)", y = "-log<sub>10</sub>(p)") + 
	 theme_bw() +
	 theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.title.x = element_markdown(size=12),axis.title.y = element_markdown(size=12),plot.title = element_text(hjust = 0.5,size=12), text = element_text(size=14)) +
	 geom_hline(yintercept=bonf_threshold, linetype="dashed", color = "orange", size=0.8) + 
	 geom_point(data = gwas_data[gwas_data$CHR==peak & -log10(gwas_data$P) > bonf_threshold & gwas_data$BP > (peak_position - 500000) & gwas_data$BP < (peak_position + 500000),], aes(x= bp_cum/1000000, y = -log10(P)),colour="orange") +
	 ggtitle(sp) 
	
	#Add butterflies images
	images = system("readlink -f /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/image_Melinaea/*", intern = TRUE)	
	my_image <- readPNG(images[1], native = TRUE)
	my_image2 <- readPNG(images[2], native = TRUE)
	print(head(gwas_data))
	
	# Set images position in the manhattan plot
	LEFT = (min(gwas_data[order(gwas_data$P)[1:3],6]) - min(gwas_data$bp_cum))/(max(gwas_data$bp_cum) - min(gwas_data$bp_cum))
	RIGHT = (max(gwas_data[order(gwas_data$P)[1:3],6]) - min(gwas_data$bp_cum))/(max(gwas_data$bp_cum) - min(gwas_data$bp_cum))
	print("Images imported")

	# Combine manhattan plot and images
	ggp_image <- p + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +                # Combine plot & image
  inset_element(p = my_image,
                left = LEFT - 0.08,
                bottom = 0.87,
                right = LEFT - 0.05,
                top = 0.90) +
  inset_element(p = my_image2,
                left = RIGHT + 0.1,
                bottom = 0.87,
                right = RIGHT + 0.13,
                top = 0.90)

	png(file=paste(sp,"_Cortex_GWAS.png",sep=""),width=600,height=400,type="cairo")
	plot(ggp_image)
	dev.off()
	
	#########################
        # GWAS zoom around peak #
        #########################
	print(paste("./gtf_extractor.sh", sp, peak_position-100000, peak_position+100000, peak, sep = " "))
	system(paste("./gtf_extractor.sh", sp, peak_position-100000, peak_position+100000, peak, sep = " "))
	print("gtf_extractor.sh ran")
	region = read.table("Gene_position.txt")
	system("rm Gene_position.txt")
	colnames(region) = c("Species","Gene","Feature","From","To","orientation","Start","End")
	start_plot = min(unlist(region[,c(7,8)]))
	end_plot = max(unlist(region[,c(7,8)]))
	print(paste(start_plot, peak_position, end_plot, peak, sep = "\t"))
	gwas_zoom = gwas_data[gwas_data$CHR==peak & gwas_data$BP > (start_plot) & gwas_data$BP < (end_plot),]
	
	p_zoom <- ggplot(gwas_zoom, aes(x = BP/1000, y = -log10(P))) +
        geom_point(color="gray") +
        scale_x_continuous() +
        labs(x = "Genome position (kb)", y = "-log<sub>10</sub>(p)") +
        theme_bw() +
        theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
	theme(axis.title.x=element_blank()) +
	theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(),plot.margin=unit(c(0.8,1,-1.2,1), "cm")) +
        geom_hline(yintercept=bonf_threshold, linetype="dashed", color = "orange", size=0.8)

	############ 
        # Syntheny #
        ############ 
	#region = region[,-c(2)]
	#region$New_start = round((region$Start + region$Stop)/2 - (region$Size/2),digits = 0)
	#region$New_end = round((region$Start + region$Stop)/2 + (region$Size/2),digits = 0)
	region_genes = region[region$Feature=="gene",]
	dummies <- make_alignment_dummies(region_genes,aes(xmin = Start, xmax = End, y = Species, id = Gene, forward=orientation),on = region_genes[1,2])
	
	# Set colors
	my_colors = rainbow(length(unique(region_genes$Gene)))
	
	# Plot syntheny based on genes only (no features plotted)
	p_syntheny <- ggplot(region_genes, aes(xmin = Start , xmax = End,y = Species, fill = Gene, forward = orientation)) +
	geom_gene_arrow()  +
	scale_color_manual(values = c("black")) + scale_fill_manual(values = my_colors) + guides(fill = guide_legend(nrow = 2)) +
	theme(legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank()) + 
	theme(panel.border = element_blank()) +
	theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"),panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="black"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,160,0), plot.margin=unit(c(1,1,-0.5,1), "cm")) 
	print("p_syntheny ready")

	# Plot syntheny based on gene and features
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
	scale_fill_manual(values = c("grey","white")) +  guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2)) +
	theme(legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank()) +
        theme(panel.border = element_blank()) +
        theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"),panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="black"),legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,0,160,0), plot.margin=unit(c(1,1,-0.5,1), "cm"))
	print("p_syntheny_features ready")

	######################################  
        # Plot syntheny and GWAS around peak #
        ######################################
	png(file=paste(sp,"_Zoom_GWAS.png",sep=""),width=600,height=500)
        grid.arrange(p_zoom,p_syntheny)
        dev.off()

	png(file=paste(sp,"_Zoom_GWAS_features.png",sep=""),width=600,height=500)
        grid.arrange(p_zoom,p_syntheny_features)
        dev.off()
	
	pdf(paste(sp,"_Zoom_GWAS_features.pdf",sep=""),9,7)
	grid.arrange(p_zoom,p_syntheny_features)
	dev.off()
	print("around peak pdf and png plotted")
	#tiff(file=paste(sp,"_Zoom_GWAS.tif",sep=""), res=300)
	#grid.arrange(p_zoom,p_syntheny)
        #dev.off()
	############################
        # Move to the next species #
        ############################
	list[[counter]] = p
	counter = counter + 1
}

###########################################
# Combine all the GWAS in a single figure #
###########################################
#png(file="Cortex_GWAS.png",width=600,height=1200,type="cairo")
#do.call(grid.arrange,c(list,ncol=1))
#dev.off()

###########################################
# Combine all the GWAS in a single figure #
###########################################
png(file="Optix_GWAS.png",width=600,height=700,type="cairo")
do.call(grid.arrange,c(list,ncol=1))
dev.off()
