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
library(gggenes)

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
	
	#####################################          
	# Extract scaffold larger than 2 Mb #
	#####################################
	list2 = list()
	a = 1
	for (i in unique(gwas_data$CHR)){
		Size = length(gwas_data[gwas_data$CHR==i,1])
		print(Size)
		if(Size > 10000){
			list2[[a]] = i
			a = a + 1
		}
	}

	Big_scaff = do.call(c,list2)
	gwas_data = gwas_data[gwas_data$CHR %in% Big_scaff,]

	########################################################################
        # Find scaffold with the peak and the associated points above the peak #
        ########################################################################
	thres_pval = 0.05/dim(gwas)[1]
	min(unlist(gwas_data$P))
	peak = gwas_data[gwas_data$P==min(unlist(gwas_data$P)),1]
	peak_position = gwas_data[gwas_data$P==min(unlist(gwas_data$P)),3]
	#peak = names(sort(table(gwas_data[gwas_data$P < thres_pval,1]),decreasing=T)[1])

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
	
	png(file=paste(sp,"_Cortex_GWAS.png",sep=""),width=600,height=400,type="cairo")
	plot(p)
	dev.off()
	
	#########################
        # GWAS zoom around peak #
        #########################
	region = read.table("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/Figure_2/Data/Genes_positions.txt")
	region = region[region$V1==sp,]
	colnames(region) = c("Species","Scaffold","Gene","Start","Stop","Size")
	start_plot = min(unlist(region[,c(4,5)]))
	end_plot = max(unlist(region[,c(4,5)]))
	gwas_zoom = gwas_data[gwas_data$CHR==peak & gwas_data$BP > (start_plot-25000) & gwas_data$BP < (end_plot+25000),]

	p_zoom <- ggplot(gwas_zoom, aes(x = bp_cum/1000000, y = -log10(P))) +
        geom_point(color="gray") +
        scale_x_continuous() +
        labs(x = "Genome position (Mb)", y = "-log<sub>10</sub>(p)") +
        theme_bw() +
        theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank(),plot.margin=unit(c(0.8,1,-1.2,1), "cm")) +
        geom_hline(yintercept=bonf_threshold, linetype="dashed", color = "orange", size=0.8)

	############ 
        # Syntheny #
        ############ 
	region = region[,-c(2)]
	region$New_start = round((region$Start + region$Stop)/2 - (region$Size/2),digits = 0)
	region$New_end = round((region$Start + region$Stop)/2 + (region$Size/2),digits = 0)
	dummies <- make_alignment_dummies(region,aes(xmin = New_start, xmax = New_end, y = Species, id = Gene),on = "parn")
	
	p_syntheny <- ggplot(posi, aes(xmin = New_start , xmax = New_end,y = Species, fill = Gene)) +
	geom_gene_arrow()  +
	scale_color_manual(values = c("black")) + scale_fill_manual(values = c("gold","azure3","coral")) + guides(color = "none")+
	theme(legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank()) + 
	theme(panel.border = element_blank()) +
	theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"),panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="black"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,200,0), plot.margin=unit(c(1,1,-0.5,1), "cm")) 
	
	######################################  
        # Plot syntheny and GWAS around peak #
        ######################################
	png(file=paste(sp,"_Zoom_GWAS.png",sep=""),width=600,height=400,type="cairo")
        grid.arrange(p_zoom,p_syntheny)
        dev.off()

	############################
        # Move to the next species #
        ############################
	list[[counter]] = p
	counter = counter + 1
}

###########################################
# Combine all the GWAS in a single figure #
###########################################
png(file="Cortex_GWAS.png",width=600,height=1200,type="cairo")
do.call(grid.arrange,c(list,ncol=1))
dev.off()
