library(dplyr)
library(readr)
library(tidyr)
library(magrittr)
library(stringr)
library(parallel)


# Read input arguments from Snakemake
input.file <- snakemake@input$complete
output.file <- snakemake@output$gene_list
cutoff <- snakemake@params$cutoff


all.samples <- read_tsv(input.file)

all.samples %>%
    # filter(grepl("^[cmnpt][LH0-9].*",sample)) %>%
    # filter(!grepl("IAA",sample)) %>%
    filter(grepl("gene", exon_id)) %>% 
    mutate(lt=length(unique(sample))) %>%
    group_by(gene_id) %>%
    summarise(mean_T_exon=sum(T_exon)/lt[1],mean_T_intron=sum(T_intron)/lt[1]) ->
    foo

all.samples %>%
    filter(gene_id %in% foo$gene_id[foo$mean_T_exon > cutoff | foo$mean_T_intron > cutoff]) %>%
    group_by(gene_id) %>%
    summarise(lt=length(unique(sample))) %>%
    ungroup %>%
    filter(lt==max(lt)) ->
    gene.list.final

write_csv(gene.list.final, path=output.file)
