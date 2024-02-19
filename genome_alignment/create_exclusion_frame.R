library(dplyr)
library(modelr)
library(readr)
library(tidyr)
library(purrr)
library(magrittr)
library(stringr)
library(cowplot)
library(superslam)
library(parallel)

read.multiple.files <- function(path,read.function,pattern,col.types=NULL,no.cores=1){
    all.files <- list.files(path=path,pattern=pattern,full.names=TRUE)
    names.all.files <- list.files(path=path,pattern=pattern,full.names=FALSE)
    all.data <- mclapply(all.files[],function(x){yy=read_tsv(x,col_names=TRUE,col_types=col.types,skip=1);colnames(yy)[ncol(yy)]="count";yy[yy$count>0,]},mc.cores=no.cores)
    all.data <- bind_rows(all.data,.id="sample")
    all.data %>%
        mutate(sample=names.all.files[as.numeric(sample)]) ->
        all.data
}



count.data.exon <- read.multiple.files(path="/extra/meisig/temp/slam_seq/intron_exon_count_comparison/",read_tsv,pattern="count.exon.*",col.types=c(.default = col_character(),Geneid=col_character(),Chr=col_character(),Start=col_double(),End=col_double()),no.cores=6)
count.data.exon %>%
    select(sample,Geneid,Length,count) ->
#    filter(count>0) ->
    count.data.exon
count.data.exon %>%
    mutate(sample=str_replace(sample,"count.","")) %>%
    mutate(sample=str_replace(sample,".txt","")) %>%
    separate(sample,into=c("type","sample","fraction"),sep="\\.") ->
    count.data.exon
count.data.exon %>%
    separate(Geneid,into=c("gene","foo"),sep="\\.") %>%
    select(-foo) ->
    count.data.exon


## count.data.intron <- read.multiple.files(path="/extra/meisig/temp/slam_seq/intron_exon_count_comparison/",read_tsv,pattern="count.intron.*",col.types=c(.default = col_character(),Geneid=col_character(),Chr=col_character(),Start=col_double(),End=col_double()),no.cores=6)
## count.data.intron %>%
##     select(sample,Geneid,Length,count) ->
## #    filter(count>0) ->
##     count.data.intron
## count.data.intron %>%
##     mutate(sample=str_replace(sample,"count.","")) %>%
##     mutate(sample=str_replace(sample,".txt","")) %>%
##     separate(sample,into=c("type","sample","fraction"),sep="\\.") ->
##     count.data.intron
## count.data.intron %>%
##     separate(Geneid,into=c("gene","foo"),sep="\\.") %>%
##     select(-foo) ->
##     count.data.intron



## count.data <- bind_rows(count.data.exon,count.data.intron)

## count.data %>%
##     group_by(type,sample) %>%
##     summarise(sum.counts=sum(count)) %>%
##     spread(type,sum.counts) %>%
##     mutate(ratio=intron/exon) ->
##     foo

## foo %>%
##     ungroup %>%
##     arrange(ratio) ->
##     foo

## foo$sample=factor(foo$sample,levels=foo$sample,ordered=T)

## foo=bind_cols(foo,convert.sample(select(foo,sample)))


## foo %>%
##     filter(!grepl("ch_",sample)) %>%
##     select(-sample) %>%
##     rename(sample=sample1)%>%
##     mutate(rep=as.numeric(rep)) %>%
## #    left_join(difference.to.nuc) %>%
##     group_by(rep,time) %>%
##     mutate(lt=sum(sample=="n")) %>%
##     ungroup %>%
##     filter(lt>0) %>%
##     filter(time>0) %>%
## #    filter(!(sample=="c" & rep==1 & time==60)) %>%
##     group_by(rep,time) %>%
##     mutate(ratio=ratio/ratio[sample=="n"]) %>%
##     ungroup %>%
##     filter(sample!="n") %>%
##     mutate(to.exclude=ratio>0.5) %>%
##     select(sample,rep,time,to.exclude) ->
##     exclusion.table


count.data.exon %>%
    select(sample,fraction) %>%
    distinct ->
    foo

foo=bind_cols(foo,convert.sample(select(foo,sample)))

foo %>%
    filter(!grepl("ch_",sample)) %>%
    filter(!grepl("IAA",sample)) %>%
    ungroup %>%
    mutate(to.exclude=FALSE) %>%
    mutate(to.exclude=ifelse(time<100|rep %in% c(5,6)|fraction!="nuc",to.exclude,TRUE)) %>%
    #intron/exon ratio
    mutate(to.exclude=ifelse(fraction!="cyto"|rep %in% c(3,4),to.exclude,TRUE)) %>%
    #intron/exon ratio
    mutate(to.exclude=ifelse(sample1!="pH"|rep %in% c(3,4),to.exclude,TRUE)) %>%
    #intron/exon ratio
    mutate(to.exclude=ifelse(sample1!="pL",to.exclude,TRUE)) %>%
    mutate(to.exclude=ifelse(sample1!="pc",to.exclude,TRUE)) %>%
    #outlier
    mutate(to.exclude=ifelse(!(sample1=="c"&time==60&rep==1),to.exclude,TRUE)) ->
    foo


foo %>%
    select(sample1,rep,time,to.exclude) %>%
    mutate(rep=as.numeric(rep)) %>%
    rename(sample=sample1) ->
    exclusion.table

write_csv(exclusion.table,path="exclusion_table_nucleus_contamination.csv")

foo %>%
    filter(to.exclude==FALSE) %>%
    select(sample) %>%
    write_csv(path="kept_samples_fitting_johannes.csv")
