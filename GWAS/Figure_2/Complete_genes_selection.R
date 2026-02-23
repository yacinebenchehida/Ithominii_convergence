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
Input = read.table("Gene_position.txt")
colnames(Input) = c("Gene","Feature","From","To","Orientation")
Input$Species = args
Input = Input[,c(6,1,2,3,4,5)]
#Input = Input %>% distinct()
output_name = "Gene_position.txt"
  
#########################################################
# Check if a gene is complete and drop incomplete genes #
#########################################################
counter = 1
list = list()
for (gn in unique(Input$Gene)){
  tmp = Input[Input$Gene==gn,]
  if(tmp$Feature[1]!="gene"){
    next
  }
  else{
    assessement = c("gene","start_codon","stop_codon") %in% unique(Input$Feature)
    assessement2 = FALSE %in% assessement
    if(assessement2 == TRUE){
      next
    }
    else if(assessement2 == FALSE){
      if(table(tmp[tmp$Feature %in% c("stop_codon","start_codon"),3])[1]==table(tmp[tmp$Feature %in% c("stop_codon","start_codon"),3])[2]){
        list[[counter]] = tmp
        counter = counter + 1
      }
      else{
        next
      }
    }
  }
}

Input = do.call(rbind,list)

####################################################################
# Make sure all genes starting their genome position at position 0 #
####################################################################
Input$start = 0
Input$end = 0

for (i in unique(Input$Gene)){
  Input[Input$Gene==i,7] = Input[Input$Gene==i & Input$Feature=="gene" ,4]
  Input[Input$Gene==i,8] = Input[Input$Gene==i & Input$Feature=="gene" ,5]
}

##############
# Set colors #
##############
my_colors = rainbow(length(unique(Input$Gene)))

###############
# orientation #
###############
Input[Input$Orientation=="+",6] = TRUE
Input[Input$Orientation=="-",6] = FALSE

#######################################
# Remove hard to visualize annotation #
#######################################
Input = Input[Input$Feature!="stop_codon",]
Input = Input[Input$Feature!="start_codon",]

#####################
# Write final table #
#####################
write.table(x = Input,file = output_name,quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)
