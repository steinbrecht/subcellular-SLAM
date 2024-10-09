library(dplyr)
library(readr)
library(tidyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
sample=args[1]
norm.data.dir=args[2]
tmp.dir=args[3]
out.dir=args[4]



print("Loading files")

file <- paste("counts_T2C.",sample,".tsv",sep="")

data <- read_tsv(paste(tmp.dir,"/",file,sep=""),col_types=c(.default = col_integer())) %>%
	  distinct

head(data)


# data %>%
#     filter(exon_Ts>0|intron_Ts>0) %>%
#     filter(grepl("ENSMUSG",gene_id)) ->
#     data


data %>%
    filter(exon_Ts>0|intron_Ts>0) %>%
    separate(gene_id,c("ID","foo"),sep="\\.") ->
    data



pconv.file <- paste("pconv_",sample,".csv",sep="")

pconv.data <- read_csv(paste(norm.data.dir,"/",pconv.file,sep=""))

print("Normalizing")
head(pconv.data)


data %>%
    mutate(sample=sample) %>%
    mutate(p_conv=pconv.data$p.conv[1]) %>%
    # separate(ID,c("gene","sub.id"),sep="_") %>%
    # select(sample,gene,sub.id,p_conv,T_exon=exon_Ts,T2C_exon=exon_T2Cs,T_intron=intron_Ts,T2C_intron=intron_T2Cs) %>%
    select(sample,ID,p_conv,T_exon=exon_Ts,T2C_exon=exon_T2Cs,T_intron=intron_Ts,T2C_intron=intron_T2Cs) %>%
    mutate(T2C_over_T_exon=(T2C_exon/T_exon)/p_conv[1],
           T2C_over_T_intron=(T2C_intron/T_intron)/p_conv[1]) ->
    data



out.file=paste("/counts_T2C_normalized.",sample,".tsv",sep="")

write_tsv(data,path=paste(out.dir,out.file,sep=""))


