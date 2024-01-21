library(ggplot2)
library(ggpubr)
library(pafr)

dat = commandArgs(trailingOnly=TRUE)
input = read.table(dat[1],header = TRUE,fill = TRUE)

ali_mars_mech <- read_paf(dat[1])
ali_mars_mech <-subset(ali_mars_mech, mapq > 7)
query = unique(ali_mars_mech$qname)
target = unique(ali_mars_mech$tname)

pdf(paste("plot_syntheny_",dat[2],"_",dat[3],".pdf",sep=""))
plot_synteny(ali_mars_mech, t_chrom=target, q_chrom=query,centre = TRUE) +
  theme_bw()
dev.off()
