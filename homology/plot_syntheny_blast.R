##################
# load libraries #
##################
library(ggplot2)
library(ggpubr)
library(pafr)

################
# Prepare data #
################
dat = commandArgs(trailingOnly=TRUE)
input = read.table(dat[1],header = TRUE,fill = TRUE)

input$querylen = dat[2]
input$subjectlen = dat[2]
input$direction = "+"
input = input[,c(seq(1,6),seq(9,12),seq(15,17))]
input = input[,c(1,2,3,11,7,8,13,4,12,9,10)]
copy= input
for (i in 1:dim(input)[1]){
  input[i,5] = input[i,5] + input[i,1] - 1
  input[i,6] = input[i,6] + input[i,1] - 1
}


#input$bad_blast = FALSE
#for (i in 2:dim(input)[1]){
#  tryCatch({
#  diff = input[i,10] - input[i-1,10]
#  print(abs(diff))
#  if(abs(diff) > 5000){
#    input[i,12] = TRUE
#    }
#  }, error=function(e){})
#}
#
#
#counter = 0
#for (i in 1:dim(input)[1]){
#  tryCatch({
#    if(input[i,12] == TRUE){
#      input[i,12] = paste("TRUE",counter,sep="")
#      counter = counter + 1
#    }
#    else{
#      counter = 0
#    }
# }, error=function(e){})
#}
#
#input = input[input$bad_blast != "TRUE0",]
#

##################################################################################
# Round value to the closest thousands so the homology blocks have the same size #
##################################################################################
input$querystart=round(input$querystart, -3)
input$queryend = input$querystart + 1000
input$subjectstart=round(input$subjectstart, -3)
input$subjectend = input$subjectstart + 1000

input$bad_blast = FALSE
for (i in 2:dim(input)[1]){
  tryCatch({
    diff = input[i,10] - input[i,5]
    print(abs(diff))
    if(abs(diff) > 5000){
      input[i,12] = TRUE
    }
  }, error=function(e){})
}

input = input[input$bad_blast != "TRUE",]
input = input[-c(1),-c(1,2)]

write.table(x = input, file = "minimap_plot.txt",quote = FALSE, row.names = FALSE,col.names = FALSE,sep="\t")


##########################################
# Plot homology using the r package pafr #
##########################################
ali_mars_mech <- read_paf("minimap_plot.txt")
query = unique(ali_mars_mech$qname)
target = unique(ali_mars_mech$tname)

pdf(paste("plot_syntheny_blast_",dat[3],"_",dat[4],".pdf",sep=""))
plot_synteny(ali_mars_mech, t_chrom=target, q_chrom=query,centre = TRUE) +
  theme_bw() +
  geom_rect(aes(xmin = 50000-43469, xmax = 50000 - 43266, ymin = 0.8 , ymax = 0.95),
            alpha = 1/2,
            fill = "red") +
  geom_rect(aes(xmin = 24750, xmax = 25250, ymin = 0.8 , ymax = 0.95),
            alpha = 1/2,
            fill = "blue") +
  geom_rect(aes(xmin = 50000-(15224528-15188153), xmax = 50000-(15224528-15187353), ymin = 2.05 , ymax = 2.20),
            alpha = 1/2,
            fill = "red") +
  geom_rect(aes(xmin = 24750, xmax = 25250, ymin = 2.05 , ymax = 2.20),
            alpha = 1/2,
            fill = "blue") 

dev.off()
