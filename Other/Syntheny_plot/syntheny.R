##################
# Load libraries #
##################
library("gggenes")
library("cowplot")
library(ggplot2)

#############
# Read data #
############# 
args = commandArgs(trailingOnly=TRUE)
Results = args[1]
Numata_genes = read.table(paste(Results,"/results.txt",sep=""))
colnames(Numata_genes) = c("Species","Gene","start","end")


#################### 
# Get genes length #
####################
genes = length(unique(Numata_genes$Gene))

####################################################################
# Make sure all genes starting their genome position at position 0 #
####################################################################
for (i in unique(Numata_genes$Species)){
  print(i)
  print(max(Numata_genes[Numata_genes$Species==i,3]))
  minimum = min(Numata_genes[Numata_genes$Species==i,3])
  print(minimum)
  Numata_genes[Numata_genes$Species==i,3] = Numata_genes[Numata_genes$Species==i,3] - minimum
  Numata_genes[Numata_genes$Species==i,4] = Numata_genes[Numata_genes$Species==i,4] - minimum
}

##########################
# Make dummies alignment #
##########################
dummies <- make_alignment_dummies(
  Numata_genes,
  aes(xmin = start, xmax = end, y = Species, id = Gene),
  on = args[2]
)

#########
# Plots #
#########
# Each species with different scale
p1 <- ggplot(Numata_genes, aes(xmin = start/1000 , xmax = end/1000, y = Species, fill = Gene)) +
  geom_gene_arrow() + scale_fill_manual(values = c(rainbow(genes))) + theme_genes() +
  facet_wrap(Species~ ., ncol = 1,scale = "free") + 
  theme(legend.position="bottom",axis.title.y=element_blank()) + 
  xlab("Relative position (kb)")

# All species with same scale
p2 <- ggplot(Numata_genes, aes(xmin = start/1000 , xmax = end/1000, y = Species, fill = Gene)) +
  geom_gene_arrow() + scale_fill_manual(values = c(rainbow(genes))) + theme_genes() +
  facet_wrap(Species~ ., ncol = 1,scale = "free_y") + 
  theme(legend.position="bottom",axis.title.y=element_blank()) + 
  xlab("Relative position (kb)")

#######################
# Save results as pdf #
#######################
pdf(paste(Results,"/syntheny_different_scale.pdf",sep=""),12,10)
p1
dev.off()

pdf(paste(Results,"/syntheny_same_scale.pdf",sep=""),12,10)
p2
dev.off()
