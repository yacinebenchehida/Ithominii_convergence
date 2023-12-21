# Packages
install.packages("stringr")
library("stringr")

# Set working directory
setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/All_Melinaea/Local_D_Stat/")

# Import data
Geno = read.table("around_cortex.genotypes.txt")
Species = read.table("Samples.txt")

# Define taxa to analyse
O = "lilis_imitata"
P1 = "mothone_mothone"
P2 = "menophilus_ssp"
P3 = "mothone_messenina"
taxa = c(P1,P2,P3,O)

# Fonction to extract calculate allelic frequencies for the 4 taxa analysed
SNP_GENOTYPE <- function(position,genotype,species){
  list=list()
  a = 1
  for(Tax in taxa){
    if(Tax==O){
      SNP = genotype[position,]
      colnames(SNP) = c("Position",species[,2])
      SNP = SNP[, (colnames(SNP) %in% Tax)]
      SNP = gsub(pattern = "[\\|\\/]",replacement = "",x = SNP)
      SNP = SNP[SNP!=".."]
      SNP = paste(SNP,collapse='')
      size = nchar(SNP)
      if(size < 4){
        break
      }
      else{
        frq_0 = str_count(SNP,"0")/size
        frq_1 = str_count(SNP,"1")/size
        list[[a]] = c(Tax,frq_0,frq_1,genotype[position,c(1)])
        a = a + 1
      }
    }
    else{
  SNP = genotype[position,]
  colnames(SNP) = c("Position",species[,2])
  SNP = SNP[, (colnames(SNP) %in% Tax)]
  SNP = gsub(pattern = "[\\|\\/]",replacement = "",x = SNP)
  SNP = SNP[SNP!=".."]
  SNP = paste(SNP,collapse='')
  size = nchar(SNP)
  if(size < 4){
    break
  }
  else{
    frq_0 = str_count(SNP,"0")/size
    frq_1 = str_count(SNP,"1")/size
    list[[a]] = c(Tax,frq_0,frq_1,genotype[position,c(1)])
    a = a + 1
      }
    }
  }
  print(list)
  SNP_genotype = as.data.frame(do.call(rbind,list))
  colnames(SNP_genotype)=c("Taxa","frq_0","frq_1","Position")
  return(SNP_genotype)
}  

# Fonction to test if a sites is monomorphic among the taxa studied
FIXED_SITES <- function(freq_data){
  if((freq_data[1,2]==1) && (freq_data[2,2]==1) && (freq_data[3,2]==1) && (freq_data[4,2]==1)){
    results = TRUE
    return(results)
  }
  else if((freq_data[1,3]==1) && (freq_data[2,3]==1) && (freq_data[3,3]==1) && (freq_data[4,3]==1)){
    results = TRUE
    return(results)
}
  else{
    results = FALSE
    return(results)
  }
}


# Compute the ABBA BABA test for each single available SNP return the ones that are fixed in the wanted disposition (or at a very high frequency)
list=list()
counter = 1

for (i in 1:dim(Geno)[1]){
  SNP_freq = SNP_GENOTYPE(i,Geno,Species)
  if(dim(SNP_freq)[1]==4){
  SNP_freq =  SNP_freq[1] |> bind_cols(sapply(SNP_freq[2:3], as.numeric)) |>  bind_cols(SNP_freq[4])
  if(apply(SNP_freq, 2, function(x) any(is.na(x)))[2]==TRUE){
     next
  }
  if(FIXED_SITES(SNP_freq)==TRUE){
    next
  }
  else if((SNP_freq[4,2]==1) || (SNP_freq[4,3]==1)){
    derived = which(SNP_freq[4,]==0)
    Pi1= SNP_freq[1,derived]
    Pi2= SNP_freq[2,derived]
    Pi3= SNP_freq[3,derived]
    Pi4= SNP_freq[4,derived]
    
    ABBA = (1-Pi1) * Pi2 * Pi3 * (1-Pi4)
    BABA = Pi1 * (1-Pi2 ) * Pi3 * (1-Pi4 )
    D_pat = (ABBA-BABA)/(ABBA+BABA)
    if(BABA!=0 && D_pat > -1){
      list[[counter]] = c(as.numeric(SNP_freq[4,4])-1100000,ABBA,BABA,D_pat)
      counter = counter + 1
     }
   }
 }
}  

patterson_D_data = data.frame(do.call(rbind,list))
colnames(patterson_D_data) = c("Position","ABBA","BABA","D")
patterson_D_data


#  Calculate number of abba/baba sites by windows of 100 kb
list=list()
counter = 1
window = 100
start = 1
end = 100
Chromosome_size =  400000

while(end < Chromosome_size){
  mid_position = (end + start)/2 - 0.5
  n_ABBA = sum(patterson_D_data[(patterson_D_data$Position >= start) & (patterson_D_data$Position < end),2])
  n_BABA = sum(patterson_D_data[(patterson_D_data$Position >= start) & (patterson_D_data$Position < end),3])
  D = mean(patterson_D_data[(patterson_D_data$Position >= start) & (patterson_D_data$Position < end),4])
  list[[counter]] = c(start+1100000,end+1100000, mid_position+1100000, n_ABBA,n_BABA,D)
  start = start + window
  end = end + window
  counter = counter + 1
}

mid_position = (Chromosome_size + start)/2 - 0.5
n_ABBA = sum(patterson_D_data[(patterson_D_data$Position >= start) & (patterson_D_data$Position < end),2])
n_BABA = sum(patterson_D_data[(patterson_D_data$Position >= start) & (patterson_D_data$Position < end),3])
D = mean(patterson_D_data[(patterson_D_data$Position >= start) & (patterson_D_data$Position < end),4])
list[[counter]] = c(start+1100000,Chromosome_size+1100000,mid_position+1100000, n_ABBA,n_BABA,D)
final_data = data.frame(do.call(rbind,list))
colnames(final_data) = c("Start","End","Mid_position","n_ABBA","n_BABA","D")
final_data

# Plotting
library(tidyr)
final_data = gather(final_data, Statistics, values, n_ABBA:D, factor_key=TRUE)

p <- ggplot(final_data, aes(x=Mid_position, y=values,color=values)) + 
  geom_point() +
  scale_colour_gradient(low = c("blue"),high=c("blue")) +
  theme_bw() + xlab("Genome position (bp)") + facet_grid(Statistics~ ., scales='free') + 
  ylab(paste("value in ", window," SNPs",sep="")) +  theme(axis.title.x = element_text(size = 25, face="bold"), axis.title.y = element_text(size = 22),axis.text.y=element_text(size=14,face="bold")) + 
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.major = element_blank()) 
plot(p)


nameP1 = gsub(pattern = "(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)", replacement = "\\1_\\2", paste(taxa,collapse = "_"))
nameP2 = gsub(pattern = "(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)", replacement = "\\3_\\4", paste(taxa,collapse = "_"))
nameP3 = gsub(pattern = "(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)", replacement = "\\5_\\6", paste(taxa,collapse = "_"))
nameO = gsub(pattern = "(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)_(\\w)\\w+_(\\w+)", replacement = "\\7_\\8", paste(taxa,collapse = "_"))
title = paste("P1:",nameP1," P2:", nameP2, " P3:", nameP3, "O:", nameO)

p <- ggplot(final_data[final_data$Statistics!="D",], aes(x=Mid_position, y=values,color=Statistics)) + 
  geom_point() +
  theme_bw() + xlab("Genome position (bp)") +
  ylab(paste("value in ", window," SNPs",sep="")) +  theme(axis.title.x = element_text(size = 25, face="bold"), axis.title.y = element_text(size = 22),axis.text.y=element_text(size=14,face="bold")) + 
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.major = element_blank()) + ggtitle(title)
plot(p)
