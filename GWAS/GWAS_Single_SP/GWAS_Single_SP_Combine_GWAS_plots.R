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
	print(head(gwas))
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
		if(Size > 5000){
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
	peak = names(sort(table(gwas_data[gwas_data$P < thres_pval,1]),decreasing=T)[1])
	print(head(gwas_data))

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
	geom_point(data = gwas_data[gwas_data$CHR==peak & -log10(gwas_data$P) > bonf_threshold,], aes(x= bp_cum/1000000, y = -log10(P)),colour="orange") +
	ggtitle(sp) 
	
	list[[counter]] = p
	counter = counter + 1
}

###########################################
# Combine all the GWAS in a single figure #
###########################################
png(file="Cortex_GWAS.png",width=600,height=1200,type="cairo")
do.call(grid.arrange,c(list,ncol=1))
dev.off()
