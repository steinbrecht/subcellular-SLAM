library(dplyr)
library(readr)
library(tidyr)
library(stringr)


args = commandArgs(trailingOnly=TRUE)
sample=args[1]
norm.data.dir=args[2]

print("Loading files")

exp.threshold <- 1000

out.file=paste("/gene_data_normalized.",sample,".csv",sep="")
gene.data <- read_csv(paste(norm.data.dir,out.file,sep=""))

gene.data %>%
    mutate(is.exonic=grepl("ENSMUSE",sub.id)) %>%
    filter(!grepl("[0-9][0-9][0-9]E",sub.id)) %>%
    filter((is.exonic & T.exon > exp.threshold)|(!is.exonic & T.intron > exp.threshold)) ->
    filtered.gene.data



out.file=paste("/gene_data_normalized_filtered.",sample,".csv",sep="")
write_csv(filtered.gene.data,path=paste(norm.data.dir,out.file,sep=""))
