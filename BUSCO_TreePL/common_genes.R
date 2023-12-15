dat = commandArgs(trailingOnly=TRUE)

busco_genes = read.table(dat,fill=TRUE)
common_genes = Reduce(intersect, busco_genes)

write.table(file = "common_genes.txt",x = common_genes,quote=FALSE,row.names=FALSE, col.names=FALSE)
