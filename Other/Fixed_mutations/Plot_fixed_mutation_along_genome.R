#!/usr/bin/env Rscript


############################
# Load necessary libraries #
############################
library(tidyr)
library(ggplot2)

##########################################
# Read argument provided externally in R #
#########################################
Argument = commandArgs(trailingOnly=TRUE)

#############################################################################################################################
# Transform external argument into the variables of interests (i.e. the data frame to read, the species and the window size #
#############################################################################################################################
dt = read.table(Argument[1], header = TRUE)
sp = Argument[2]
Window_size = as.numeric(Argument[3])/1000

############################################################################
# Linear the genome position (otherwise it starts at zero for each contig) #
############################################################################
dt$genome_position = 0
for (i in 1:dim(dt)[1]){
  if(i==1){
    dt[i,6] = dt[i,3]
  }
  else{
    dt[i,6] = dt[i-1,6] + (dt[i,3]-dt[i,2]) + 1
  }
}

##########################################################################################################################
# Transform input data frame from wide to long format (combining together the two types of mutations (fixed and het/hom) #
########################################################################################################################## 
a = gather(dt, Type_mutations, value, Fixed:Het_Hom, factor_key=TRUE)

#####################################################
# Plot the two types of mutations across the genome #
#####################################################
pdf(paste(Argument[4],"/Fixed_het_hom_mutations_",sp,".pdf",sep=""),7,3.5)
ggplot(a, aes(x=genome_position/1000000, y=value,fill=Type_mutations)) + 
  geom_point(size = 1.5,colour="black",shape = 21) + 
  scale_fill_manual(name = "Fixed\nmutations",values = c("purple","orange")) +
  theme_bw() + xlab("Genome position (Mb)") +
  ylab(paste("Mutations per ", Window_size,"kb windows",sep="")) + facet_grid(.~Type_mutations) +
  theme(legend.box.background = element_rect(colour = "black",size = 1)) + ggtitle(gsub("_"," ",sp)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
