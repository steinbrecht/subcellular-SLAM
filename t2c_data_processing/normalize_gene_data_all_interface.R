library(dplyr)
library(readr)
library(tidyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
sample=args[1]
norm.data.dir=args[2]
muttable.dir=args[3]
out.dir=args[4]



print("Loading files")

freq.file <- paste("gene_T2C.",sample,".tsv",sep="")

gene.data <- read_tsv(paste(muttable.dir,"/",freq.file,sep=""),col_types=c(.default = col_integer())) %>%
	  distinct

head(gene.data)


gene.data %>%
    filter(exon_Ts>0|intron_Ts>0) %>%
    filter(grepl("ENSMUSG",gene_id)) ->
    gene.data


gene.data %>%
    filter(exon_Ts>0|intron_Ts>0) %>%
    separate(gene_id,c("gene","foo"),sep="\\.") ->
    gene.data



pconv.file <- paste("pconv_",sample,".csv",sep="")

pconv.data <- read_csv(paste(norm.data.dir,"/",pconv.file,sep=""))

print("Normalizing")
head(pconv.data)


gene.data %>%
    mutate(sample=sample) %>%
    mutate(p.conv=pconv.data$p.conv[1]) %>%
    separate(gene,c("gene","sub.id"),sep="_") %>%
    select(sample,gene,sub.id,p.conv,T.exon=exon_Ts,T2C.exon=exon_T2Cs,T.intron=intron_Ts,T2C.intron=intron_T2Cs) %>%
    mutate(t2c.over.t.exon=(T2C.exon/T.exon)/p.conv[1],t2c.over.t.intron=(T2C.intron/T.intron)/p.conv[1]) ->
    gene.data



out.file=paste("/gene_data_normalized.",sample,".csv",sep="")

write_csv(gene.data,path=paste(out.dir,out.file,sep=""))


