library(GenomicRanges)
library(seqQTL)
library(tidyverse)
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

remove.exon.entry <- function(gene.irange,exon.irange,gene.grange,exon.grange,split.hits,x){
    rv <- GenomicRanges::setdiff(gene.irange[as.numeric(x),],exon.irange[split.hits[[x]]])
    rv <- GRanges(seqnames=rep(seqnames(gene.grange)[1],length(rv)),ranges=rv)
    rv$gene=gene.grange[as.numeric(x),]$gene
    rv$exon=exon.grange[split.hits[[x]]]$exon
    rv$type="exon"
    rv$parent=paste(gene.grange[as.numeric(x),]$gene,"_T",sep="")
    rv
}

retain.exon.info <- function(dta){
    gene.start <- min(dta[dta$X3=="exon",]$X4)
    gene.end <- max(dta[dta$X3=="exon",]$X5)
    gene=paste(dta$gene.id[1],"_gene",sep="")
    gene.tbl=tibble(type="gene",parent=NA,start=gene.start,end=gene.end,gene=gene,exon=NA)
    transcript.tbl=tibble(type="transcript",parent=gene,start=gene.start,end=gene.end,gene=gene,exon=NA)
    exon.tbl <- tibble(type="exon",parent=paste(gene,"T",sep="_"),start=dta[dta$X3=="exon",]$X4,end=dta[dta$X3=="exon",]$X5,gene=gene,exon=dta[dta$X3=="exon",]$exon.id,count=dta[dta$X3=="exon",]$count.avg)
    bind_rows(gene.tbl,transcript.tbl,exon.tbl)
}

rdc.helper <- function(to.reduce,reducer,ignore.strand=F){
    hits <- findOverlaps(to.reduce, reducer,ignore.strand=ignore.strand)
    grl <- extractList(reducer, as(hits, "List"))
    psd <- pintersect(to.reduce, grl,ignore.strand=ignore.strand)
    multiples <- unlist(lapply(psd,length))
                                        #    new.metadata=unlist(lapply(1:length(psd),function(x){rep(to.reduce$gene[x],multiples[x])}))
    new.metadata=bind_rows(lapply(1:length(psd),function(x){(as.data.frame(mcols(to.reduce)[rep(x,multiples[[x]]),]))}))
    psd=unlist(psd)
    mcols(psd)=new.metadata
                                        #uniquify because
                                        #....---....... reducer
                                        #...------..... w1
                                        #..------...... w2
                                        #intersection w1,reducer = intersection w2,reducer
                                        #....---....... result
    
    psd=psd[!duplicated(ranges(psd)),]
                                        #    psd <- lapply(1:length(psd),function(x){if (length(psd[[x]])>0){rv=psd[[x]]$gene=to.reduce[x,]$gene}else{rv=psd[[x]]};rv})
    psd
}


slice.gene.new <- function(info.list,window.width,window.shift){
    dta=info.list[["original"]]
    intron.restriction=info.list[["intron.windows"]]
    gene.start <- min(dta[dta$X3=="exon",]$X4)
    gene.end <- max(dta[dta$X3=="exon",]$X5)
    if (gene.end-gene.start<window.width){
        borders=cbind(gene.start,gene.end)
    }else{
        no.shifts <- floor((gene.end-gene.start-window.width)/window.shift)
        borders <- cbind(gene.start+window.shift*0:no.shifts,gene.start+window.width+window.shift*0:no.shifts)
        borders[nrow(borders),2] <- gene.end
    }
    gene.names=paste(dta$gene.id[1],as.character(sprintf("%05d",1:nrow(borders))),sep="_")
    gene.ranges.object <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=as.numeric(borders[,1]),end=as.numeric(borders[,2])),gene=gene.names,type="gene")
    ab <- GRanges(seqnames="foo",ranges=IRanges(start=intron.restriction$start,end=intron.restriction$end))
    gene.ranges.object <- rdc.helper(gene.ranges.object,ab)
    gene.ranges.object=gene.ranges.object[width(gene.ranges.object)>window.width/4,]
#    gene.ranges.object <- unlist(gene.ranges.object)
                                        #    transcripts <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=as.numeric(borders[,1]),end=as.numeric(borders[,2])),gene=gene.names,parent=gene.names,type="transcript")
    if (length(gene.ranges.object)>0){
    gene.names=paste(dta$gene.id[1],as.character(sprintf("%05d",1:length(gene.ranges.object))),sep="_")
    gene.ranges.object$gene=gene.names
        transcripts <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=start(gene.ranges.object),end=end(gene.ranges.object)),gene=gene.ranges.object$gene,parent=gene.ranges.object$gene,type="transcript")
    }else{
        transcripts <- GRanges(seqnames=NULL,strand=NULL,ranges=IRanges(start=NULL,end=NULL),gene=NULL,parent=NULL,type=NULL)
        gene.ranges.object <- GRanges(seqnames=NULL,strand=NULL,ranges=IRanges(start=NULL,end=NULL),gene=NULL,parent=NULL,type=NULL)
    }
#    exon.ranges.object <- GRanges(seqnames="foo",strand=dta$X7[1],IRanges(start=dta[dta$X3=="exon",]$X4,end=dta[dta$X3=="exon",]$X5),exon=dta[dta$X3=="exon",]$exon.id)
#    exon.irange <- IRanges(start=start(exon.ranges.object),end=end(exon.ranges.object))
#    gene.irange <- IRanges(start=start(gene.ranges.object),end=end(gene.ranges.object))
#    exon.irange <- GenomicRanges::reduce(exon.irange)
#    exon.ranges.object$exon=paste(as.character(sprintf("%05d",1:length(exon.ranges.object))),"E",sep="")
#    ovl <- as_tibble(findOverlaps(exon.irange,gene.irange,type="any"))
#    split.hits <- split(ovl$queryHits,ovl$subjectHits)
#    new.gene.grange <- lapply(names(split.hits),remove.exon.entry,gene.irange=gene.irange,exon.irange=exon.irange,gene.grange=gene.ranges.object,exon.grange=exon.ranges.object,split.hits=split.hits)
#    exon.list <- lapply(names(split.hits),generate.exon.entry.new,gene.irange=gene.irange,exon.irange=exon.irange,gene.grange=gene.ranges.object,exon.grange=exon.ranges.object,split.hits=split.hits)
#    exon.range <- do.call("c",exon.list)
#    return.value <- as_tibble(c(gene.ranges.object,transcripts,exon.range))[,c("type","parent","start","end","gene","exon")]
    if (length(transcripts)>0){
        return.value <- as_tibble(c(gene.ranges.object,transcripts))[,c("type","parent","start","end","gene")]
    }else{
        return.value <- tibble(type=character(),parent=character(),start=integer(),end=integer(),gene=character())
}
return.value
}


assemble.X9 <- function(x){
    if (x$X3=="gene"){
        id.string=paste("ID=",x$gene,";","gene_id=",x$gene,sep="")
    }
    if (x$X3=="transcript"){
        id.string=paste("ID=",x$gene,"_T",";","gene_id=",x$gene,";","Parent=",x$parent,sep="")
    }
    if (x$X3=="exon"){
                                        #        id.string=paste("ID=",x$exon,";","gene_id=",x$gene,";","Parent=",x$parent,sep="")
        true.gene.id=strsplit(x$gene,split="_")[[1]][1]
        id.string=paste("ID=",x$exon,";","exon_id=",true.gene.id,"_",x$exon,";","gene_id=",x$gene,";","Parent=",x$parent,";","Count=",sprintf("%.2f",x$count),sep="")
    }
                                        #    tibble(X9=id.string)
    id.string
}



args = commandArgs(trailingOnly=TRUE)
count.file=args[1]
gff.file=args[2]
no.cores=as.numeric(args[3])
out.path=args[4]
count.file="/extra/meisig/temp/slam_seq/results_exon_new/count_data/count_all_samples.csv"
gff.file="/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3"
no.cores=4
counts <- read_csv(count.file)


counts %>%
#    filter(fraction %in% c("nuc","cyto")) %>%
    filter(fraction == "nuc") %>%
    group_by(fraction,gene,exon,Chr,Start,End,Strand,Length) %>%
    summarise(count.avg=mean(count)) ->
    counts

counts %>%
    ungroup %>%
    select(gene,exon,count.avg) %>%
    rename(gene.id=gene,exon.id=exon) ->
    counts


#count.threshold <- 20

## counts %>%
##     filter(fraction == "nuc") %>%
##     group_by(fraction,gene,exon,Chr,Start,End,Strand,Length) %>%
##     summarise(count.avg=mean(count)) %>%
##     filter(count.avg>count.threshold) ->
##     counts.filtered


## counts.filtered %>%
##     group_by(gene,exon) %>%
##     mutate(lt=length(fraction)) %>%
## #    filter(lt==2) %>%
##     dplyr::select(gene,exon,Chr,Start,End,Strand,Length) %>%
##     distinct ->
##     counts.filtered


#gff <- read_tsv("/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3",col_names=F,guess_max=1e4,quote="",skip=7,n_max=1e4)
gff <- read_tsv(gff.file,col_names=F,guess_max=1e4,quote="",skip=7)

chr.names <- paste("chr",c(as.character(1:19),"X","Y","M"),sep="")

gff[,] %>%
    filter(X3 %in% c("gene","exon")) ->
    gff

ids <- mclapply(gff$X9,extract.ids,mc.cores=no.cores)

ids <- bind_rows(ids)
gff <- bind_cols(gff,ids)

gff %>%
    left_join(counts) ->
    gff

## gff[,] %>%
##     filter(is.na(exon.id)|exon.id %in% counts.filtered$exon) ->
##     gff.filtered

## gff.filtered %>%
##     group_by(gene.id) %>%
##     mutate(no.exons=length(na.omit(unique(exon.id)))) %>%
##     ungroup %>%
##     filter(no.exons>0) ->
##     gff.filtered

#write_csv(gff.filtered,"~/projektentwicklung/slam_seq/new_annotation_assembly/gene_id.gff")

## gff.filtered %>%
##     group_by(gene.id) %>%
##     group_split() ->
##     split.gff.filtered

## gff.filtered %>%
##     group_by(gene.id) %>%
##     group_keys() ->
##     split.gff.filtered.annot


gff %>%
    group_by(gene.id) %>%
    group_split() ->
    split.gff

gff %>%
    group_by(gene.id) %>%
    group_keys() ->
    split.gff.annot


#names(split.gff.filtered) <- split.gff.filtered.annot$gene.id
names(split.gff) <- split.gff.annot$gene.id


load(file="intron_mask.RData")
#split.gff.filtered <- split.gff.filtered[names(split.gff.filtered) %in% names(intron.gff)]
split.gff <- split.gff[names(split.gff) %in% names(intron.gff)]


intron.window.info <- lapply(names(split.gff),function(x){list(original=split.gff[[x]],intron.windows=intron.gff[[x]])})

#subgroup.ix <- sample(1:length(split.gff),100)
#subgroup.ix <- 1:length(split.gff.filtered)

#retain.exon.info(split.gff[[sample(1:length(split.gff),1)]])

#slice.gene.new(split.gff[[sample(1:length(split.gff),1)]],window.width=2e3,window.shift=1e2)
new.annotation.exon <- mclapply(split.gff[],retain.exon.info,mc.cores=no.cores)
#new.annotation <- mclapply(split.gff.filtered[],slice.gene.new,window.width=2e3,window.shift=1e2,mc.cores=no.cores)
#new.annotation <- mclapply(intron.window.info[],slice.gene.new,window.width=2e3,window.shift=1e2,mc.cores=no.cores)
new.annotation <- mclapply(intron.window.info[],slice.gene.new,window.width=4e3,window.shift=4e2,mc.cores=no.cores)

new.annotation <- bind_rows(bind_rows(new.annotation),bind_rows(new.annotation.exon))

new.annotation %>%
    separate(gene,into=c("gene.id","window"),sep="_",remove=F) %>%
    select(type,gene.id,start,end,gene,exon,parent,count) ->
    new.annotation

filter(gff,X3=="gene") %>%
    dplyr::select(X1,X2,X6,X7,X8,X9,gene.id) %>%
    distinct ->
    chromosome.and.strand.information


new.annotation %>%
    dplyr::select(type,gene.id,start,end,gene,exon,parent,count) %>%
    left_join(chromosome.and.strand.information,by="gene.id") %>%
    dplyr::select(X1,X2,X3=type,X4=start,X5=end,X6,X7,X8,X9,gene,exon,parent,count) ->
    new.annotation.full

assembled.X9 <- mclapply(1:nrow(new.annotation.full),function(x){assemble.X9(new.annotation.full[x,])},mc.cores=no.cores)

new.annotation.full$X9 <- unlist(assembled.X9)
new.annotation.out <- new.annotation.full[,paste("X",as.character(1:9),sep="")]
#write_tsv(new.annotation.out,path="gencode_derived_sliding_window_genes.gff",col_names=F)

new.annotation.out %>%
    group_by(X1) %>%
    group_split ->
    annotation.list

names(annotation.list) <- unlist(lapply(annotation.list,function(x){x$X1[1]}))

for (i in 1:length(annotation.list)){
    new.annotation.out <- annotation.list[[i]]
    write_tsv(new.annotation.out,path=paste(out.path,"/sliding_intron_parallel_",names(annotation.list)[i],".gff",sep=""),col_names=F)
}

## gff <- read_tsv("~/projektentwicklung/slam_seq/new_annotation_assembly/sliding_intron_all_chr.gff",col_names=F)

## extract.ids <- function(string){
##     gene.id <- str_extract(string,"ENSMUSG[0-9]*_[0-9][0-9][0-9]")
##     list(gene.id=gene.id)
## }


## gff[,] %>%
##     filter(X3=="gene") %>%
##     group_by(X1,X2,X3,X4,X5,X6,X7,X8,X9) %>%
##     do(as_tibble(extract.ids(.$X9))) %>%
##     mutate(gene=gene.id) %>%
##     ungroup %>%
##     select(gene,X4,X5,X7) ->
##     gff

## write_csv(gff,path="~/projektentwicklung/slam_seq/new_annotation_assembly/sliding_intron_reduced_all_chr.gff")

