library(GenomicRanges)
library(seqQTL)
library(tidyverse)
library(broom)
library(magrittr)

extract.ids <- function(string){
    gene.id <- str_extract(string,"ENSMUSG[0-9]*")
    exon.id <- str_extract(string,"ENSMUSE[0-9]*")
    transcript.id <- str_extract(string,"ENSMUST[0-9]*")
    list(gene.id=gene.id,exon.id=exon.id,transcript.id=transcript.id)
}


generate.exon.entry.new <- function(gene.irange,exon.irange,gene.grange,exon.grange,split.hits,x){
    rv <- GenomicRanges::intersect(gene.irange[as.numeric(x),],exon.irange[split.hits[[x]]])
    rv <- GRanges(seqnames=rep(seqnames(gene.grange)[1],length(rv)),ranges=rv)
    rv$gene=gene.grange[as.numeric(x),]$gene
    rv$exon=exon.grange[split.hits[[x]]]$exon
    rv$type="exon"
    rv$parent=paste(gene.grange[as.numeric(x),]$gene,"_T",sep="")
    rv
}

slice.gene.new <- function(dta,window.width,window.shift){
    gene.start <- min(dta[dta$X3=="exon",]$X4)
    gene.end <- max(dta[dta$X3=="exon",]$X5)
    if (gene.end-gene.start<window.width){
        borders=cbind(gene.start,gene.end)
    }else{
        no.shifts <- floor((gene.end-gene.start-window.width)/window.shift)
        borders <- cbind(gene.start+window.shift*0:no.shifts,gene.start+window.width+window.shift*0:no.shifts)
        borders[nrow(borders),2] <- gene.end
    }
    gene.names=paste(dta$gene.id[1],as.character(sprintf("%03d",1:nrow(borders))),sep="_")
    gene.ranges.object <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=as.numeric(borders[,1]),end=as.numeric(borders[,2])),gene=gene.names,type="gene")
    transcripts <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=as.numeric(borders[,1]),end=as.numeric(borders[,2])),gene=gene.names,parent=gene.names,type="transcript")

    exon.ranges.object <- GRanges(seqnames="foo",strand=dta$X7[1],IRanges(start=dta[dta$X3=="exon",]$X4,end=dta[dta$X3=="exon",]$X5),exon=dta[dta$X3=="exon",]$exon.id)
    exon.irange <- IRanges(start=start(exon.ranges.object),end=end(exon.ranges.object))
    gene.irange <- IRanges(start=start(gene.ranges.object),end=end(gene.ranges.object))
    exon.irange <- GenomicRanges::reduce(exon.irange)
    exon.ranges.object$exon=as.character(sprintf("%03d",1:length(exon.ranges.object)))
    ovl <- as_tibble(findOverlaps(exon.irange,gene.irange,type="any"))
    split.hits <- split(ovl$queryHits,ovl$subjectHits)
    exon.list <- lapply(names(split.hits),generate.exon.entry.new,gene.irange=gene.irange,exon.irange=exon.irange,gene.grange=gene.ranges.object,exon.grange=exon.ranges.object,split.hits=split.hits)
    exon.range <- do.call("c",exon.list)
    return.value <- as_tibble(c(gene.ranges.object,transcripts,exon.range))[,c("type","parent","start","end","gene","exon")]
}


assemble.X9 <- function(x){
    if (x$X3=="gene"){
        id.string=paste("ID=",x$gene,";","gene_id=",x$gene,sep="")
    }
    if (x$X3=="transcript"){
        id.string=paste("ID=",x$gene,"_T",";","gene_id=",x$gene,";","Parent=",x$parent,sep="")
    }
    if (x$X3=="exon"){
        id.string=paste("ID=",x$exon,";","gene_id=",x$gene,";","Parent=",x$parent,sep="")
    }
                                        #    tibble(X9=id.string)
    id.string
}


counts <- read_csv("/extra/meisig/temp/slam_seq/results_exon_new/count_data/count_all_samples.csv")


count.threshold <- 20

counts %>%
    filter(fraction %in% c("nuc","cyto")) %>%
    group_by(fraction,gene,exon,Chr,Start,End,Strand,Length) %>%
    summarise(count.avg=mean(count)) %>%
    filter(count.avg>count.threshold) ->
    counts.filtered

counts.filtered %>%
    group_by(gene,exon) %>%
    mutate(lt=length(fraction)) %>%
    filter(lt==2) %>%
    dplyr::select(gene,exon,Chr,Start,End,Strand,Length) %>%
    distinct ->
    counts.filtered


#gff <- read_tsv("/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3",col_names=F,guess_max=1e4,quote="",skip=7,n_max=1e4)
gff <- read_tsv("/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3",col_names=F,guess_max=1e4,quote="",skip=7)

chr.names <- paste("chr",c(as.character(1:19),"X","Y","M"),sep="")

gff[,] %>%
    filter(X3 %in% c("gene","exon")) %>%
    group_by(X1,X2,X3,X4,X5,X6,X7,X8,X9) %>%
    do(as_tibble(extract.ids(.$X9))) %>%
    filter(is.na(exon.id)|exon.id %in% counts.filtered$exon) ->
    gff

gff %>%
    group_by(gene.id) %>%
    mutate(no.exons=length(na.omit(unique(exon.id)))) %>%
    ungroup %>%
    filter(no.exons>1) ->
    gff

write_csv(gff,"~/projektentwicklung/slam_seq/new_annotation_assembly/gene_id.gff")

gff %>%
    group_by(gene.id) %>%
    do(slice.gene.new(.,window.width=2e3,window.shift=1e2)) ->
    new.annotation

filter(gff,X3=="gene") %>%
    dplyr::select(X1,X2,X6,X7,X8,X9,gene.id) %>%
    distinct ->
    chromosome.and.strand.information


new.annotation %>%
    dplyr::select(type,gene.id,start,end,gene,exon,parent) %>%
    left_join(chromosome.and.strand.information,by="gene.id") %>%
    dplyr::select(X1,X2,X3=type,X4=start,X5=end,X6,X7,X8,X9,gene,exon,parent) ->
    new.annotation.full


new.annotation.full[,] %>%
    group_by(gene,exon,X1,X2,X3,X4,X5,X6,X7,X8) %>%
    do(X9=assemble.X9(.)) %>%
    tidy(X9) ->
    new.annotation.full.small

colnames(new.annotation.full.small)[ncol(new.annotation.full.small)] <- "X9"

new.annotation.out <- new.annotation.full.small[,paste("X",as.character(1:9),sep="")]
#write_tsv(new.annotation.out,path="gencode_derived_sliding_window_genes.gff",col_names=F)


new.annotation.out %>%
    group_by(X1) %>%
    do(ll=(.)) %$%
    setNames(ll,X1) ->
    annotation.list

for (i in 1:length(annotation.list)){
    new.annotation.out <- annotation.list[[i]]
    write_tsv(new.annotation.out,path=paste("sliding_intron_",names(annotation.list)[i],".gff",sep=""),col_names=F)
}

gff <- read_tsv("~/projektentwicklung/slam_seq/new_annotation_assembly/sliding_intron_all_chr.gff",col_names=F)

extract.ids <- function(string){
    gene.id <- str_extract(string,"ENSMUSG[0-9]*_[0-9][0-9][0-9]")
    list(gene.id=gene.id)
}


gff[,] %>%
    filter(X3=="gene") %>%
    group_by(X1,X2,X3,X4,X5,X6,X7,X8,X9) %>%
    do(as_tibble(extract.ids(.$X9))) %>%
    mutate(gene=gene.id) %>%
    ungroup %>%
    select(gene,X4,X5,X7) ->
    gff

write_csv(gff,path="~/projektentwicklung/slam_seq/new_annotation_assembly/sliding_intron_reduced_all_chr.gff")

