# Packages 
library(ggplot2)
library(dplyr)
library(reshape2)
library("gplots")
library(wesanderson)
library(gtools)
library(ComplexHeatmap)

# read data 
args = commandArgs(trailingOnly=TRUE)
Path_Data = args[1]
data = read.table(Path_Data, colClasses = c(rep("factor",5)))

# Prepare data in the right format
data = data[,-c(1)]
colnames(data) = c("Sample_name","Phenotype","Genotype","Species","CHR","SNP")
data$Genotype = gsub(pattern = "/", x = data$Genotype,replacement = "",perl = TRUE)
data$Genotype = gsub(pattern = "\\|", x = data$Genotype,replacement = "",perl = TRUE)

# Read color data
Path_Data_Color = args[2]
color = read.table(Path_Data_Color,header = TRUE,comment.char='@',encoding="UTF-8")
couleur = color[unique(data$Species) %in% color$Species,]
data$Species = gsub(pattern = "(\\w+)_(\\w+)", x = data$Species,replacement = "\\2",perl = TRUE)
data$CHR = gsub(pattern = "SUPER_", x = data$CHR,replacement = "Chr ",perl = TRUE)
data = data[order(data$Species),]
couleur$Species = gsub(pattern = "(\\w+)_(\\w+)", x = color$Species,replacement = "\\2",perl = TRUE)

# Extract information about the phenotype for the heatmap annotation
a = table(data$Phenotype)/length(unique(data$SNP))
a = a[a!=0]
list=list()
list2=list()
counter = 1
color = c("orange","sienna","darkgoldenrod1")

for (i in 1:length(a)){
  list[[counter]]= rep(names(a[i]),a[i])
  bizbizfaitlamouche_leretour = color[counter]
  names(bizbizfaitlamouche_leretour)= c(names(a[i]))
  list2[[counter]]= bizbizfaitlamouche_leretour
  counter = counter + 1
}


# Extract information about the species for the heatmap annotation
list3=list()
list4=list()
b = table(data$Species)/length(unique(data$SNP))
pheno1 = unique(data[which(data$Phenotype=="2"),4])
pheno2 = unique(data[which(data$Phenotype=="3"),4])
pheno_intermediate = unique(data[which(data$Phenotype=="2.5"),4])
b = c(b[pheno1],b[pheno_intermediate],b[pheno2])
counter = 1

for (i in 1:length(b)){
  list3[[counter]]= rep(names(b[i]),b[i])
  bizbizfaitlamouche_leretour = couleur[couleur$Species == names(b[i]),2]
  names(bizbizfaitlamouche_leretour)= c(names(b[i]))
  list4[[counter]]= bizbizfaitlamouche_leretour
  counter = counter + 1
}



# Extract information about the CHROMOSOME for the heatmap annotation
chr = table(data$CHR)/length(unique(data$Sample_name))
list7=list()
list8=list()
counter = 1
color = rep("black",length(chr))

for (i in 1:length(chr)){
  list7[[counter]]= rep(names(chr[i]),chr[i])
  bizbizfaitlamouche_leretour = color[counter]
  names(bizbizfaitlamouche_leretour)= c(names(chr[i]))
  list8[[counter]]= bizbizfaitlamouche_leretour
  counter = counter + 1
}

# Chromosome information for plotting
chromosome_tab = table(data$CHR)/length(unique(data$Sample_name))
list_chr = list()
counter = 1
for (i in 1:length(chromosome_tab)){
  list_chr[[counter]] = rep(names(chromosome_tab[i]),chromosome_tab[i])
  counter = counter + 1
}

chromo_plot = do.call(c,list_chr)
chromo_plot = mixedsort(chromo_plot)


# Extract information about the SNPS for the heatmap annotation
pos = table(unique(data$SNP))
list5=list()
list6=list()
counter = 1
color = rainbow(length(unique(data$SNP)))

for (i in 1:length(pos)){
  list5[[counter]]= rep(names(pos[i]),pos[i])
  bizbizfaitlamouche_leretour = color[counter]
  names(bizbizfaitlamouche_leretour)= c(names(pos[i]))
  list6[[counter]]= bizbizfaitlamouche_leretour
  counter = counter + 1
}

# reorder data so species fits the genotype
list_sp = list()
counter = 1
for (i in names(b)){
  list_sp[[counter]] = data[data$Species==i,]
  counter = counter + 1
}
data = as.data.frame(do.call(rbind, list_sp))



# Information plotted in the heatmap
list_snps = list()
counter = 1
for (i in unique(data$SNP)){
  list_snps[[counter]] = data[data$SNP==i,3]
  counter = counter + 1
}
Input = as.matrix(do.call(cbind, list_snps))

# Set the specific annotation block used by complexheatmap
Prefix = args[3]
Number_pheno = length(levels(data$Phenotype))

if(Number_pheno==2){
	Phenotypes_names = c(Prefix, paste("No ",Prefix, sep = ""))
	ha_row = rowAnnotation(Pheno=do.call(c,list), Species = do.call(c,list3), col = list(Pheno = do.call(c,list2), Species =  do.call(c,list4)),border = TRUE,annotation_legend_param = list(Pheno = list(labels = Phenotypes_names)))
}else{
	Phenotypes_names = c(Prefix, paste("No ",Prefix, sep = ""),"Hybrid")
	ha_row = rowAnnotation(Pheno=do.call(c,list), Species = do.call(c,list3), col = list(Pheno = do.call(c,list2), Species =  do.call(c,list4)),border = TRUE,annotation_legend_param = list(Pheno = list(labels = Phenotypes_names)))
}

ha_column = HeatmapAnnotation(Chromosome = do.call(c,list7), col = list( Chromosome = do.call(c,list8)),border = TRUE)

 
# Define color code
pal <- wes_palette("Zissou1", length(unique(data$Genotype)), type = "continuous")
pal = c("white",pal[c(1,3,4)])

# Plot heatmap
pdf(paste(Prefix,"_GWAS_heatmap.pdf",sep=""),7,8)
Heatmap(Input,show_row_dend=FALSE,show_column_dend=FALSE,name = "zygosity",show_row_names = FALSE,use_raster=FALSE, row_title_gp = gpar(fontsize = 10), bottom_annotation = ha_column, left_annotation = ha_row,row_split= factor(do.call(c,list)), column_split = factor(chromo_plot,levels = mixedsort(levels(factor(chromo_plot)))), col = pal,border = TRUE)
dev.off()
