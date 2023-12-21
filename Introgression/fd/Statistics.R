library(tidyr)
library(ggplot2)

# Set working directory
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

# Import data
Prefix=args[2]
data = read.table(paste(args[1],"/Results/",Prefix,"/",args[3],sep=""),header = T)
window=as.numeric(args[4])

# Define taxa to analyse
P1 = args[5] 
P2 = args[6] 
P3 = args[7] 
O = args[8] 
taxa = c(P1,P2,P3,O)

# Prepare data
data[data$D < 0, 10] = 0
data[data$fd < 0, 10 ] = 0
data[data$fd > 1, 10 ] = 0
data[data$fdM < 0, 11 ] = NA
data = data[data$sitesUsed > window/10 ,]
final_data = gather(data, Statistics, values, ABBA:fdM, factor_key=TRUE)
head(final_data)


# Add the species analysed in the plot
nameP1 = gsub(pattern = "(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)", replacement = "\\1_\\2", paste(taxa,collapse = "_"))
nameP2 = gsub(pattern = "(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)", replacement = "\\3_\\4", paste(taxa,collapse = "_"))
nameP3 = gsub(pattern = "(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)", replacement = "\\5_\\6", paste(taxa,collapse = "_"))
nameO = gsub(pattern = "(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)", replacement = "\\7_\\8", paste(taxa,collapse = "_"))
title = paste("P1:",nameP1," P2:", nameP2, " P3:", nameP3, "O:", nameO)



# Plot data
pdf(paste(args[1],"/Results/",Prefix,"/",Prefix,"_All_stats_plot.pdf",sep=""))
p <- ggplot(final_data, aes(x=mid, y=values,color=values)) + 
  geom_point() +
  theme_bw() + xlab("Genome position (bp)") +
  ylab(paste("value in ", window," SNPs",sep="")) +  theme(axis.title.x = element_text(size = 25, face="bold"), axis.title.y = element_text(size = 22),axis.text.y=element_text(size=14,face="bold")) + 
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.major = element_blank())  + facet_grid(Statistics~ ., scales='free')
plot(p)
dev.off()

pdf(paste(args[1],"/Results/",Prefix,"/",Prefix,"_ABBA_BABA_stats_plot.pdf",sep=""))
p <- ggplot(final_data[!(final_data$Statistics %in% c("D","fdM","fd")),], aes(x=mid, y=values,color=Statistics)) + 
  geom_point() +
  theme_bw() + xlab("Genome position (bp)") +
  ylab(paste("value in ", window," SNPs",sep="")) +  theme(axis.title.x = element_text(size = 25, face="bold"), axis.title.y = element_text(size = 22),axis.text.y=element_text(size=14,face="bold")) + 
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.major = element_blank()) + ggtitle(title)
plot(p)
dev.off()
