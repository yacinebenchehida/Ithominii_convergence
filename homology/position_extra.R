dat = commandArgs(trailingOnly=TRUE)
data = read.table(dat)

scaffold = data$V1[1]
start = data$V2[2] - 150000
end = data$V2[2] + 150000

cat(paste(scaffold, start, end, sep="\t"),"\n")
