library(admixtools)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
subspecies_files <- args[1]
data = args[2]
pops <- read.table(subspecies_files)
pop1 <- pops[1,1]
pop2 <- pops[2,1]
pop3 <- pops[3,1]
pop4 <- pops[4,1]
outfile = paste(args[2],"/",pop1,"_",pop2,"_",pop3,"_",pop4,".txt", sep="")
print(outfile)

out = f4(data, pop1, pop2, pop3, pop4)

write_tsv(out, outfile)

