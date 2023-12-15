library(ape)
library(ggtree)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

#############
# Load data #
#############
Prefix = args[2]

myTree <- ape::read.tree(paste(args[1],"/",Prefix,"_Rooted_final_output.newick",sep=""))
data = read.table("/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/Twisst/8_Single_tree/Inputs/Grouping.txt")
colnames(data) = c("Ind","Species")

#############################################
# Convert individual labels to species name #
#############################################
for (i in data$Ind){
  myTree$tip.label[which(grepl(i, myTree$tip.label))] = data[data$Ind==i,2]
}

##################################################################
# Write phylogeny with species names instead of individual names #
##################################################################
write.tree(myTree, file = paste(args[1],"/",Prefix,".raxml_renamed.bestTree.newick",sep=""), append = FALSE, digits = 5, tree.names = FALSE)

##################
# Set color code #
##################
#couleur = c("chartreuse3","aquamarine1","aquamarine","forestgreen","forestgreen","black","black","blue","dodgerblue","mediumblue","sienna","coral","red","salmon","brown1","brown","darkorange","gold","orchid","darkorchid","#8F4783","purple","darkviolet","grey")
couleur = c("blue","red","forestgreen")
colors = as.data.frame(cbind(unique(sort(myTree$tip.label)),couleur))
colnames(colors) = c("Species","Color")

#################################################################
# Assign to each individual the color associated to its species #
#################################################################
data$Color = "2BeAdded"
for (i in 1:dim(data)[1]) {
  for (j in 1:dim(colors)[1]) {
    if(data[i,2]==colors[j,1]){
      data[i,3] = colors[j,2]
    }  
  }
}

##########
# ggtree #
##########
tree2 <- ggtree(myTree)
data_ggtree = as.data.frame(cbind(tree2$data$label,"2BeAdded"))

for (i in 1:dim(data_ggtree)[1]) {
  for (j in 1:dim(colors)[1]){
    tryCatch({
    if(data_ggtree[i,1]==colors[j,1]){
      data_ggtree[i,2] = colors[j,2]
    }
  },error=function(e){})
 }
}

tree2$data$col <- cut(as.numeric(tree2$data$label),
                      breaks = c(-1, 70, 90, 100),
                      labels = c("<70",">70", ">90"))

p <- ggtree(tree2$data) + geom_tiplab(size =4, color= data_ggtree[data_ggtree$V2!="2BeAdded",2]) 
p <- p + geom_point2(aes(subset=(!isTip), fill= tree2$data$col, size = tree2$data$col, shape=tree2$data$col))  + theme(legend.position="right")  + scale_fill_manual(breaks = c("<70",">70", ">90") , values = c("green","gold","red")) + 
  scale_size_manual(breaks = c("<70",">70", ">90") , values = c(1.2,2,2.5)) + scale_shape_manual(breaks = c("<70",">70", ">90") , values = c(21,21,21)) + guides(size = FALSE, shape=FALSE, fill = guide_legend(override.aes=list(shape=21))) + labs(fill="Bootstraps") 
p = p + geom_treescale(x=0.31) + ggtitle(args[2])

pdf(paste(args[1],"/",Prefix,"_tree.pdf",sep=""),25,13) 
p
dev.off()
