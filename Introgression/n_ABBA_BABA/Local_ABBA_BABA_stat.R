library("stringr")
library("ggplot2")
library("dplyr")

# Set working directory
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

# Import data
Prefix=args[2]
Geno = read.table(paste(args[1],"/Results/",Prefix,".genotypes.txt",sep=""))
Species = read.table(paste(args[1],"/Data/Samples.txt",sep=""))
cat("data loaded\n")

# Define taxa to analyse and window size
P1 = args[3] 
P2 = args[4] 
P3 = args[5] 
O = args[6] 
taxa = c(P1,P2,P3,O)
Start = as.numeric(args[7])
End = as.numeric(args[8])
cat("taxa loaded\n")

# Fonction to extract calculate allelic frequencies for the 4 taxa analysed
SNP_GENOTYPE <- function(position,genotype,species){
  list=list()
  a = 1
  for(Tax in taxa){
    if(Tax==O){ # If it's the outgroup
      SNP = genotype[position,]
      colnames(SNP) = c("Position",species[,2])
      SNP = SNP[, (colnames(SNP) %in% Tax)]
      SNP = gsub(pattern = "[\\|\\/]",replacement = "",x = SNP)
      SNP = SNP[SNP!=".."] # Remove samples that with missing data
      SNP = paste(SNP,collapse='') # Transform vector into a simple character
      size = nchar(SNP) # Count number of characters
      frq_0 = str_count(SNP,"0")/size # Frequency of 1st allele
      frq_1 = str_count(SNP,"1")/size # Frequency of the 2nd allele
      list[[a]] = c("Outgroup",frq_0,frq_1,genotype[position,c(1)]) 
      a = a + 1
    }
    else{ # If it's not the outgroup
  SNP = genotype[position,]
  colnames(SNP) = c("Position",species[,2])
  SNP = SNP[, (colnames(SNP) %in% Tax)]
  SNP = gsub(pattern = "[\\|\\/]",replacement = "",x = SNP)
  SNP = SNP[SNP!=".."]
  SNP = paste(SNP,collapse='')
  size = nchar(SNP)
  frq_0 = str_count(SNP,"0")/size
  frq_1 = str_count(SNP,"1")/size
  list[[a]] = c(Tax,frq_0,frq_1,genotype[position,c(1)])
  a = a + 1
    }
  }
  SNP_genotype = as.data.frame(do.call(rbind,list))
  colnames(SNP_genotype)=c("Taxa","frq_0","frq_1","Position")
  return(SNP_genotype)
}  

cat("function allele frequencies ready\n")

# Fonction to test if a sites is monomorphic among the taxa studied
FIXED_SITES <- function(freq_data){
  if((freq_data[1,2]==1) && (freq_data[2,2]==1) && (freq_data[3,2]==1) && (freq_data[4,2]==1)){ # If all individuals are 0 at the SNP
    results = TRUE
    return(results)
  }
  else if((freq_data[1,3]==1) && (freq_data[2,3]==1) && (freq_data[3,3]==1) && (freq_data[4,3]==1)){ # If all individuals are 1 at the SNP 
    results = TRUE
    return(results)
}
  else{
    results = FALSE
    return(results)
  }
}

cat("function fixed sites ready\n")


# Compute the ABBA BABA test for each single available SNP return the ones that are fixed in the wanted disposition (or at a very high frequency)
list=list()
counter = 1
print(Geno)

for (i in 1:dim(Geno)[1]){
  SNP_freq = SNP_GENOTYPE(i,Geno,Species)
  SNP_freq =  SNP_freq[1] |> bind_cols(sapply(SNP_freq[2:3], as.numeric)) |>  bind_cols(SNP_freq[4]) # Make sure the SNPs values are numeric not characters
  if(apply(SNP_freq, 2, function(x) any(is.na(x)))[2]==TRUE){ # If there is a nan or na value (i.e. only missing values at least for one taxa at the SNP) skip the SNP.
     next
  }
  if(FIXED_SITES(SNP_freq)==TRUE){ # if the SNP is fixed skip the SNP
    next
  }
  else if((SNP_freq[4,2]==1) || (SNP_freq[4,3]==1)){ # If the SNP is not polymorphic in the outgroup
    derived = which(SNP_freq[4,]==0) # Find the derived position. Assume that the outgroup allele is the ancestral one.
    Pi1= SNP_freq[1,derived] # Extract allele frequency of the derived allele at P1
    Pi2= SNP_freq[2,derived] # Extract allele frequency of the derived allele at P2
    Pi3= SNP_freq[3,derived] # Extract allele frequency of the derived allele at P3
    Pi4= SNP_freq[4,derived] # Extract allele frequency of the derived allele at O (which is always 0)
    
    ABBA = Pi3*(1-Pi4)*((1-Pi1)*Pi2-Pi1*(1-Pi2)) # Calculate ABBA based on allelic frequencies
    BABA = Pi3*(1-Pi4)*((1-Pi1)*Pi2+Pi1*(1-Pi2)) # Calculate BABA based on allelic frequencies
    
    D_pat = ABBA/BABA # Compute the Patterson's D as the ratio between the two
    if(BABA!=0 && D_pat > 0.80){ # if D is larger than 0.80 (meaning P2 and P3 are closer)
      list[[counter]] = c(as.numeric(SNP_freq[4,4])-Start,D_pat)
      counter = counter + 1
    }
  }
}

patterson_D_data = data.frame(do.call(rbind,list))
colnames(patterson_D_data) = c("Position","D")
patterson_D_data


#  Calculate number of abba/baba sites by windows of 100 kb
list=list()
counter = 1
window = 100
start = 1
end = 100
Chromosome_size =  End - Start

while(end < Chromosome_size){
  mid_position = (end + start)/2 - 0.5
  Sites = dim(patterson_D_data[(patterson_D_data$Position >= start) & (patterson_D_data$Position < end),])[1]/100
  list[[counter]] = c(start+Start,end+Start, mid_position+Start, Sites)
  start = start + window
  end = end + window
  counter = counter + 1
}

mid_position = (Chromosome_size + start)/2 - 0.5
Sites = dim(patterson_D_data[(patterson_D_data$Position >= start) & (patterson_D_data$Position < end),])[1]
list[[counter]] = c(start+Start,Chromosome_size+Start,mid_position+Start, Sites)
final_data = data.frame(do.call(rbind,list))
colnames(final_data) = c("Start","End","Mid_position","D")
final_data

# Plotting
p <- ggplot(final_data, aes(x=Mid_position, y=D,color=D)) + 
  geom_point() +
  scale_colour_gradient(name = "Patterson's D",low = c("blue"),high=c("blue")) +
  theme_bw() + xlab("Genome position (bp)") +
  ylab(paste("High D per 100 SNPs",sep="")) +  theme(axis.title.x = element_text(size = 25, face="bold"), axis.title.y = element_text(size = 22),axis.text.y=element_text(size=14,face="bold")) + 
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.major = element_blank()) 

pdf(paste(args[1],"/Results/",Prefix,"_plot.pdf",sep=""))
plot(p)
dev.off()
