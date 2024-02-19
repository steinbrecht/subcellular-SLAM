library(dplyr)
library(readr)
library(tidyr)
library(magrittr)
library(stringr)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
input.file=args[1]
output.file=args[2]

all.samples <- read_csv(input.file)

all.samples %>%
    filter(grepl("^[cmnpt][LH0-9].*",sample)) %>%
    filter(!grepl("IAA",sample)) %>%
    mutate(lt=length(unique(sample))) %>%
    group_by(gene) %>%
    summarise(mean.T.exon=sum(T.exon)/lt[1],mean.T.intron=sum(T.intron)/lt[1]) ->
    foo

all.samples %>%
    filter(gene %in% foo$gene[foo$mean.T.exon>1e3]) %>%
    group_by(gene) %>%
    summarise(lt=length(unique(sample))) %>%
    ungroup %>%
    filter(lt==max(lt)) ->
    gene.list.final

write_csv(gene.list.final,path=output.file)
