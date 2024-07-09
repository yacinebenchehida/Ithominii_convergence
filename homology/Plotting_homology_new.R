#################
# Read raw data #
#################
library(dplyr)
library(ggplot2)
library(viridis)

dat = commandArgs(trailingOnly=TRUE)


input = read.table(dat[1],header = TRUE,fill = TRUE)
query_species = dat[2]
target_species = dat[3]
#target_centering =  input[1,5]
#subject_centering = input[1,10]

####################################
# "Linearise" the query coordinates #
####################################
for (i in 1:dim(input)[1]){
  input[i,5] = input[i,5] + input[i,1] - 1
  input[i,6] = input[i,6] + input[i,1] - 1

}

##############################################
# Keep only the most reliable alignment bits #
##############################################
input = input[input$Quality > 40,]
input = input[,-c(1,2)]


#########################################################
# Reformat the data in the geom_polygon expected format #
#########################################################
list = list()
counter = 1

for (i in 1:dim(input)[1]){
  list[[counter]] <- data.frame(
    x = c(input$queryStart[i],input$queryEnd[i],  input$SubjectEnd[i],input$SubjectStart[i]),
    y = c(1,1,3,3),
    group = paste("window_",i,sep=""), 
    order = c(1, 2, 4, 3),
    reference = c(query_species, query_species, target_species, target_species))
  counter = counter + 1
}

plotting_data <- do.call(rbind,list)

###################################
# Load the annotation information #
###################################
anno <- read.table("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Inputs/annotation.txt",sep="\t")
anno <- anno[anno$V1 %in% c(query_species,target_species),]


# Function to check the coordinate of the annotation and adjust the value. Also if the reference genome is reverse complemented it reorder the annotation so it fits the default orientation.
check_direction_and_realign_annotation <- function(anno){
  if(anno$V4[2] - anno$V4[1] > 0){
    star_position <- min(anno$V4)
    anno$V4 <- anno$V4 - star_position
    anno$V5 <- anno$V5 - star_position
    return(anno)
  }
  else{
    star_position <- min(anno$V4)
    anno$V4 <- anno$V4 - star_position
    anno$V5 <- anno$V5 - star_position
    anno$size <- anno$V5 - anno$V4
    anno$succes_size <- 0
    for (i in 1:(dim(anno)[1]-1)){
      anno$succes_size[i] <- anno$V4[i] - anno$V5[i+1]
    }
    anno$V4[1] = 0
    anno$V5[1] =  anno$size[1]
    for (i in 2:(dim(anno)[1])){
      anno$V4[i] <- anno$V5[i-1] + anno$succes_size[i-1]
      anno$V5[i] <- anno$V4[i] + anno$size[i]
    }
    return(anno[,c(1:5)])
  }
}

#########################################################
# Prepare the properly annotation file for the plotting #
##########################################################
list = list()
a = 1
for (i in c(query_species,target_species)){
  tmp <- anno[anno$V1==i,]
  tmp <- check_direction_and_realign_annotation(tmp)
  tmp <- tmp[tmp$V5< max(input$queryEnd),]
  list[[a]] = tmp
  a = a + 1
}

anno = as.data.frame(do.call(rbind,list))

########################################################################
# offset the annotation so the coordinate fit the one of the alignment #
########################################################################
anno[anno$V1==unique(anno$V1)[1],4]= anno[anno$V1==unique(anno$V1)[1],4] + input$queryStart[1] 
anno[anno$V1==unique(anno$V1)[1],5]= anno[anno$V1==unique(anno$V1)[1],5] + input$queryStart[1]
anno[anno$V1==unique(anno$V1)[2],4]= anno[anno$V1==unique(anno$V1)[2],4] + input$SubjectStart[1] 
anno[anno$V1==unique(anno$V1)[2],5]= anno[anno$V1==unique(anno$V1)[2],5] + input$SubjectStart[1]

#####################################################################################################
# Extract alignment regions in the annotation file. These bits will be colorised in bed in the plot #
#####################################################################################################
to_colorise <- anno[anno$V1==query_species,c(4,5)]

list = list()
counter = 1
for (i in unique(plotting_data$group)){
  tmp = plotting_data[plotting_data$group==i & plotting_data$y==1,]
  for (j in 1:dim(to_colorise)[1]){
  if(any(tmp$x >=  to_colorise[j,1] & tmp$x <= to_colorise[j,2])){
          list[[counter]] = plotting_data[plotting_data$group==i,]
          counter = counter + 1
      }
    }
  }

to_colorise <- do.call(rbind,list)

############
# Plotting #
############
pdf(paste("plot_syntheny_windows_",dat[2],"_",dat[3],".pdf",sep=""))
ggplot(plotting_data, aes(x = x, y = y, group = group)) +
  geom_polygon(fill="black",alpha = 0.5) +
  theme_minimal() + 
  geom_rect(aes(xmin = 0, xmax = max(plotting_data[plotting_data$y=="1",1]) , ymin = 0.95, ymax = 0.75), colour = "black", fill = "white") + 
  scale_y_continuous(breaks = unique(plotting_data$y), labels = c("menophilus","mothone")) +
  geom_rect(data = anno[anno$V1==unique(anno$V1)[1],], aes(xmin = V4, xmax = V5, ymin = 0.948 , ymax = 0.752,  fill = V3),   inherit.aes = FALSE) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = viridis(5)) +
  geom_rect(aes(xmin = input$SubjectStart[1] , xmax = max(plotting_data[plotting_data$y=="3",1]), ymin = 3.05, ymax = 3.25), colour = "black", fill = "white") +
  geom_rect(data = anno[anno$V1==unique(anno$V1)[2],], aes(xmin = V4, xmax = V5, ymin = 3.052 , ymax = 3.248, fill = V3),inherit.aes = FALSE) +
  theme(axis.ticks.x = element_blank()) +
  geom_polygon(data=to_colorise, aes(x = x, y = y, group = group), fill="red",alpha = 0.5) +
  labs(x = NULL, y = NULL) 
dev.off()
