library(dplyr)
library(ggplot2)
library(viridis)


dat = commandArgs(trailingOnly=TRUE)

aln = read.table(dat[1],header = TRUE,fill = TRUE)
query_species = dat[3]
target_species = dat[2]
aln <- aln[aln$Identity > 60,]


###################################
# Load the annotation information #
###################################
anno <- read.table("/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/homology/Inputs/annotation.txt",sep="\t")
anno <- anno[anno$V1 %in% c(query_species,target_species),]

# Function to check the coordinate of the annotation and adjust the value. Also if the reference genome is reverse complemented it reorder the annotation so i$
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
  #  tmp <- tmp[tmp$V5< max(input$queryEnd),]
  list[[a]] = tmp
  a = a + 1
}

anno = as.data.frame(do.call(rbind,list))

############################
# Figure out genome offset #
############################
offset <- abs(floor((aln$queryStart[1]-aln$start[1]) / 1000) * 1000)

for (i in 1:dim(aln)[1]){
  aln[i,9] = aln[i,9] + aln[i,1] - 1 
  aln[i,10] = aln[i,10] + aln[i,1] - 1 
  aln[i,5] = aln[i,5] 
  aln[i,6] = aln[i,6] 
}

##########################################################################################
# Remove interval aligning to very distant non homologous part of the alternative genome #
##########################################################################################                                          
size_difference <- abs(max(anno[anno$V1==target_species,5])-max(anno[anno$V1==query_species,5]))
if(size_difference < 5000){
  size_difference = 5000
}

tolerance <-  seq(from = 5000, to = size_difference + 2000, length.out = dim(aln)[1])
list = list()
counter = 1

differential <- aln$SubjectEnd[length(aln$SubjectEnd)] - aln$queryEnd[length(aln$queryEnd)]

for (i in 1:dim(aln)[1]){
  tol <- tolerance[i]
  difference <- aln$SubjectEnd[i]-aln$queryEnd[i]
       if(abs(difference) <  tol & differential < 0 & difference < 0){
         print(paste(i,abs(difference), differential, tol))
          list[[counter]] = aln[i,]
          counter = counter + 1
        }
      else if(abs(difference) <  tol & differential > 0 & difference > 0){
      list[[counter]] = aln[i,]
      counter = counter + 1
      }
    else if(abs(difference) < 5000){
      list[[counter]] = aln[i,]
      counter = counter + 1
    }
  }

aln <- as.data.frame(do.call(rbind,list))

#########################################################
# Reformat the data in the geom_polygon expected format #
#########################################################
list = list()
counter = 1

for (i in 1:dim(aln)[1]){
  list[[counter]] <- data.frame(
    x = c(aln$queryStart[i],aln$queryEnd[i],  aln$SubjectEnd[i],aln$SubjectStart[i]),
    y = c(1,1,3,3),
    group = paste("window_",i,sep=""), 
    order = c(1, 2, 4, 3),
    reference = c(query_species, query_species, target_species, target_species))
  counter = counter + 1
}

plotting_data <-  do.call(rbind,list)

#####################################################################################################
# Extract alignment regions in the annotation file. These bits will be colorised in bed in the plot #
#####################################################################################################
to_colorise_query <- anno[anno$V1==query_species,c(3,4,5)]
to_colorise_target <- anno[anno$V1==target_species,c(3,4,5)]

list = list()
counter = 1
for (i in unique(plotting_data$group)){
  tmp1 = plotting_data[plotting_data$group==i & plotting_data$y==1,]
  tmp2 = plotting_data[plotting_data$group==i & plotting_data$y==3,]
  for (j in 1:dim(to_colorise_query)[1]){
    if(any(tmp1$x >=  to_colorise_query[j,2] & tmp1$x <= to_colorise_query[j,3])){
      if(any(tmp2$x >=  to_colorise_target[j,2] & tmp2$x <= to_colorise_target[j,3])){
      list[[counter]] = plotting_data[plotting_data$group==i,]
      counter = counter + 1
      }
    }
  }
}

to_colorise <- do.call(rbind,list)


############
# Plotting #
############
pdf(paste("plot_syntheny_windows_nucmer_",dat[2],"_",dat[3],".pdf",sep=""),12,9)
ggplot(plotting_data, aes(x = x, y = y, group = group)) +
  geom_polygon(fill="black",alpha = 0.3) +
  theme_minimal() + 
  geom_rect(aes(xmin = 0, xmax = max(anno[anno$V1==query_species,5]) , ymin = 0.95, ymax = 0.75), colour = "black", fill = "white") + 
  scale_y_continuous(breaks = unique(plotting_data$y), labels = c(query_species,target_species)) +
  geom_rect(data = anno[anno$V1==unique(anno$V1)[1],], aes(xmin = V4, xmax = V5, ymin = 0.948 , ymax = 0.752,  fill = V3),   inherit.aes = FALSE) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = viridis(6)) +
  geom_rect(aes(xmin = 0, xmax = max(anno[anno$V1==target_species,5]), ymin = 3.05, ymax = 3.25), colour = "black", fill = "white") +
  geom_rect(data = anno[anno$V1==unique(anno$V1)[2],], aes(xmin = V4, xmax = V5, ymin = 3.052 , ymax = 3.248, fill = V3),inherit.aes = FALSE) +
  theme(axis.ticks.x = element_blank()) +
  geom_polygon(data=to_colorise, aes(x = x, y = y, group = group), color="brown",fill="brown") +
  labs(x = NULL, y = NULL) 
dev.off()
