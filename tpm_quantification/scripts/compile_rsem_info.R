library(dplyr)
library(tidyr)
library(readr)
library(parallel)
library(stringr)


convert.sample <- function(x){
    x%>%
    mutate(sample=str_replace(sample,"min","_")) %>%
    mutate(sample=str_replace(sample,"ch_","")) %>%
    separate(sample,c("sample","rep"),sep="_") %>%
    mutate(time=str_extract(sample,"(\\d)+")) %>%
    mutate(sample=str_replace(sample,"[0-9].","")) %>%
    mutate(sample=str_replace(sample,"[0-9]","")) %>%
    mutate(time=as.numeric(time))
}


read.multiple.files <- function(path,read.function,pattern,col.types=NULL,no.cores=1){
    all.files <- list.files(path=path,pattern=pattern,full.names=TRUE)
    names.all.files <- list.files(path=path,pattern=pattern,full.names=FALSE)
    all.data <- mclapply(all.files,read.function,col_names=TRUE,col_types=col.types,mc.cores=no.cores)
    all.data <- bind_rows(all.data,.id="sample")
    all.data %>%
        mutate(sample=names.all.files[as.numeric(sample)]) ->
        all.data
}

args = commandArgs(trailingOnly=TRUE)
rsem.dir=args[1]


options(readr.num_columns=0)
all.rsem <- read.multiple.files(rsem.dir,read.function="read_tsv",pattern="rsem_t.*.isoforms.results")

all.rsem %>%
    filter(!grepl("IAA",sample)) %>%
    separate(sample,into=c("foo","sample","bar"),sep="_") %>%
    select(-foo,-bar) ->
    all.rsem

all.rsem %>%
    select(sample,transcript_id,gene_id,TPM) %>%
    spread(sample,TPM) ->
    rsem.spread

all.rsem %>%
    filter(sample!="t60min5") ->
    all.rsem

all.rsem %>%
    group_by(transcript_id,gene_id) %>%
    summarise(TPM=mean(TPM)) ->
    rsem.mean

rsem.mean %>%
    filter(TPM>0) %>%
    separate(gene_id,into=c("gene","foo"),sep="\\.") %>%
    separate(transcript_id,into=c("transcript","foo"),sep="\\.") %>%
    select(-foo) %>%
    group_by(gene) %>%
    top_n(1,TPM) %>%
    group_by(gene) %>%
    filter(transcript==transcript[1]) ->
    rsem.top

write_csv(rsem.top,path="rsem_top_expressed_transcripts_total.csv")


all.rsem <- read.multiple.files(rsem.dir,read.function="read_tsv",pattern="*.isoforms.results")

all.rsem %>%
    filter(!grepl("IAA",sample)) %>%
    filter(!grepl("ch_",sample)) %>%
    separate(sample,into=c("foo","sample","bar"),sep="_") %>%
    select(-foo,-bar) ->
    all.rsem

all.rsem %>%
    select(sample,transcript_id,gene_id,TPM) %>%
    convert.sample %>%
    dplyr::select(sample,time,rep,gene_id,transcript_id,TPM) ->
    all.rsem

all.rsem %>%
    filter(sample %in% c("n","c","m","t")) %>%
    filter(TPM>0) %>%
    separate(rep,c("rep","foo","bar"),sep="\\.") %>%
    separate(gene_id,into=c("gene","foo"),sep="\\.") %>%
    separate(transcript_id,into=c("transcript","foo"),sep="\\.") %>%
    select(-foo,-bar) %>%
    filter(transcript %in% rsem.top$transcript) %>%
    group_by(gene) %>%
    filter(transcript==transcript[1]) %>%
    spread(sample,TPM) %>%
    mutate(ratio=c/n) %>%
    write_csv("tpm_normalized_fraction_counts_not_averaged.csv")

all.rsem %>%
    group_by(sample,transcript_id,gene_id) %>%
    summarise(TPM=mean(TPM)) ->
    rsem.mean

rsem.total=read_csv("rsem_top_expressed_transcripts_total.csv")

rsem.mean %>%
    filter(TPM>0) %>%
    separate(gene_id,into=c("gene","foo"),sep="\\.") %>%
    separate(transcript_id,into=c("transcript","foo"),sep="\\.") %>%
    select(-foo) %>%
    filter(transcript %in% rsem.total$transcript) %>%
    group_by(gene) %>%
    filter(transcript==transcript[1]) ->
    rsem.top


rsem.top %>%
#    select(-transcript)%>%
    spread(sample,TPM) ->
    out.frame


write_csv(out.frame,path="tpm_normalized_fraction_counts.csv")

