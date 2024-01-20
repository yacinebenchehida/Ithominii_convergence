library(ggplot2)
library(ggpubr)
library(pafr)

dat = commandArgs(trailingOnly=TRUE)
input = read.table(dat[1],header = TRUE,fill = TRUE)


for (i in 1:dim(input)[1]){
  input[i,4] =  input[i,9] 
  input[i,5] = input[i,5] + input[i,1] - 1
  input[i,6] = input[i,6] + input[i,1] - 1
}


input = input[input$Quality > 7,]
input = input[,-c(1,2)]

write.table(x = input, file = "minimap_plot.txt",quote = FALSE, row.names = FALSE,col.names = FALSE,sep="\t")

ali_mars_mech <- read_paf("minimap_plot.txt")
query = unique(ali_mars_mech$qname)
target = unique(ali_mars_mech$tname)

pdf(paste("plot_syntheny_",dat[2],"_",dat[3],".pdf",sep=""))
plot_synteny(ali_mars_mech, t_chrom=target, q_chrom=query,centre = TRUE) +
  theme_bw()
dev.off()
