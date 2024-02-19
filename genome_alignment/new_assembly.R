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

create.new.annot <- function(id.list){
    type=c("gene","transcript","exon")
    id.string=paste(id.list[["gene.id"]],id.list[["exon.number"]],c("G","T","E"),sep="_")
    parent.string=c("",paste("Parent=",id.string[1],sep=""),paste("Parent=",id.string[2],sep=""))
    id.string=paste("ID=",id.string,sep="")
    tibble(type=type,id=id.string,parent=parent.string,gene.id=rep(paste("gene_id=",paste(id.list[["gene.id"]],id.list[["exon.number"]],sep="_"),sep=""),3))
}

create.new.annot <- function(id.list){
    no.exons=length(id.list$exon.ids)
    id.list$exon.ids <- sprintf("%04d",id.list$exon.ids)
    type=c("gene","transcript",rep("exon",no.exons))
    id.string.gene.transcript=c(id.list[["gene.id"]],paste(id.list[["gene.id"]],"T",sep="_"))
    id.string.all=c(id.string.gene.transcript,paste(id.list[["gene.id"]],id.list[["exon.ids"]],sep="_"))
    parent.string.gene.transcript=c("",paste("Parent=",id.list[["gene.id"]],sep=""))
    parent.string.all=c(parent.string.gene.transcript,rep(paste("Parent=",id.string.gene.transcript[2],sep=""),no.exons))
    id.string.all=paste("ID=",id.string.all,sep="")
    list(id.string.all,parent.string.all)
    tibble(type=type,id=id.string.all,parent=parent.string.all,gene.id=rep(id.list[["gene.id"]],length(id.string.all)))
    #gene.id=rep(paste("gene_id=",paste(id.list[["gene.id"]],id.list[["exon.number"]],sep="_"),sep=""),3))
}



#gff <- read_tsv("/extra/meisig/temp/slam_seq/reference/Mus_musculus.GRCm38.92_flat.gtf",col_names=F,guess_max=1e6,quote="")
## gff <- read_tsv("/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3",col_names=F,guess_max=1e6,quote="",skip=7)

## gff <- gff %>%
##     filter(X3=="exon")

## gff[,] %>%
##     group_by(X1,X2,X3,X4,X5,X6,X7,X8,X9) %>%
##     do(as_tibble(extract.ids(.$X9))) ->
##     gff


## cc <- colnames(gff)
## cc[grepl(".id",cc)] <- c("gene","exon","transcript")
## colnames(gff) <- cc

## gff <- ungroup(gff)

## GRanges(seqnames=counts.filtered$Chr,ranges=IRanges(start=counts.filtered$Start,end=counts.filtered$End),strand=counts.filtered$Strand,gene=counts.filtered$gene,exon=counts.filtered$exon) ->
##     exons.filtered

## ss <- split(exons.filtered, ~gene)
## gene.names <- names(ss)
## aa <- reduce(split(exons.filtered, ~gene))
## new.names.genes <- rep(names(aa),elementNROWS(aa))
## new.names.exons <- unlist(lapply(elementNROWS(aa),function(x){1:x}))
## aa <- unlist(aa)
## aa$gene=new.names.genes
## aa$exon=new.names.exons
## #aa$gene=new.names
## #GR2gff(regions=aa, "blub.gff", feature.type = "experimental_feature", src = "GenomicRanges", score = ".", phase = ".")


## create.new.annot(list(gene.id=aa$gene[1],exon.number=aa$exon[1])) %>% as.data.frame


## create.new.annot(list(gene.id=aa$gene[1],exon.ids=aa$exon[1:3])) %>% 
##     unite("X9",id,parent,gene.id,sep=";") %>%
##     mutate(X9=str_replace(X9,";;",";"))

## as_tibble(aa)[,] %>%
##     group_by(seqnames,gene,start,end,width,strand) %>%
##     do(create.new.annot(list(gene.id=.$gene[1],exon.ids=.$exon))) %>%
##     unite("X9",id,parent,gene.id,sep=";") %>%
##     ungroup %>%
##     mutate(X9=str_replace(X9,";;",";")) %>%
##     mutate(X2="manual_assembly",X6="NA",X8="NA") %>%
##     dplyr::select(seqnames,X2,type,start,end,X6,strand,X8,X9) ->
##     gff.out

## colnames(gff.out) <- paste("X",1:9,sep="")

## counts <- read_csv("/extra/meisig/temp/slam_seq/results_exon_new/count_data/count_all_samples.csv")


## counts[,] %>%
##     ungroup %>%
##     left_join(gff,by=c("gene","exon")) %>%
##     select(fraction,sample,time,rep,gene,exon,transcript,Start,End,Strand,Length,count) ->
##     counts.annotated

## count.threshold <- 20

## counts %>%
##     filter(fraction %in% c("nuc","cyto")) %>%
##     group_by(fraction,gene,exon,Chr,Start,End,Strand,Length) %>%
##     summarise(count.avg=mean(count)) %>%
##     filter(count.avg>count.threshold) ->
##     counts.filtered

## counts.filtered %>%
##     group_by(gene,exon) %>%
##     mutate(lt=length(fraction)) %>%
##     filter(lt==2) %>%
##     dplyr::select(gene,exon,Chr,Start,End,Strand,Length) %>%
##     distinct ->
##     counts.filtered

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
gff <- read_tsv("/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3",col_names=F,guess_max=1e4,quote="",skip=7) %>% filter(X1=="chrM")

chr.names <- paste("chr",c(as.character(1:19),"X","Y"),sep="")

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

gff %>%
    group_by(X1) %>%
    mutate(ss=sample(1:length(gene.id),length(gene.id)))

gff %>%
    group_by(X1) %>%
    mutate(gg=paste(sample(unique(gene.id))[1:3],collapse="_")) %>%
    group_by(gene.id) %>%
    mutate(ll=grepl(gene.id[1],gg[1])) %>%
    filter(ll==TRUE|gene.id=="ENSMUSG00000012396")->
    gff



           

generate.exon.entry <- function(a1,a2,split.hits,x){
    rv <- GenomicRanges::intersect(a1[as.numeric(x),],a2[split.hits[[x]]])
    rv$gene=a1[as.numeric(x),]$gene
    rv$exon=a2[split.hits[[x]]]$exon
    rv$type="exon"
    rv$parent=paste(a1[as.numeric(x),]$gene,"_T",sep="")
    rv
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
    exons.start <- dta[dta$X3=="exon",]$X4
    gene.start <- min(exons.start)
    exon.end <- max(dta[dta$X3=="exon",]$X5)
    gene.end <- exon.end
    dta[dta$X3=="gene",]$X5 <- gene.end
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
#    exon.range$exon=paste(exon.range$gene,exon.range$exon,sep="_")
    return.value <- as_tibble(c(gene.ranges.object,transcripts,exon.range))[,c("type","parent","start","end","gene","exon")]
}

slice.gene <- function(dta,window.width,window.shift){
#print(dta$gene.id[1])
    ## if (dta$X7=="-"){
    ##     foo <- dta$X4
    ##     bar <- dta$X5
    ##     dta$X4 <- bar
    ##     dta$X5 <- foo
    ## }
#    gene.start <- dta[dta$X3=="gene",]$X4
    exons.start <- dta[dta$X3=="exon",]$X4
    gene.start <- min(exons.start)
    exon.end <- max(dta[dta$X3=="exon",]$X5)
    gene.end <- exon.end
    dta[dta$X3=="gene",]$X5 <- gene.end
    if (gene.end-gene.start<window.width){
        borders=cbind(gene.start,gene.end)
    }else{
        no.shifts <- floor((gene.end-gene.start-window.width)/window.shift)
        borders <- cbind(gene.start+window.shift*0:no.shifts,gene.start+window.width+window.shift*0:no.shifts)
        borders[nrow(borders),2] <- gene.end
    }
    #print(borders)
#    print(borders[,2]-borders[,1])
    gene.names=paste(dta$gene.id[1],as.character(sprintf("%03d",1:nrow(borders))),sep="_")
    gene.ranges.object <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=as.numeric(borders[,1]),end=as.numeric(borders[,2])),gene=gene.names,type="gene")
    transcripts <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=as.numeric(borders[,1]),end=as.numeric(borders[,2])),gene=gene.names,parent=gene.names,type="transcript")

    exon.ranges.object <- GRanges(seqnames="foo",strand=dta$X7[1],IRanges(start=dta[dta$X3=="exon",]$X4,end=dta[dta$X3=="exon",]$X5),exon=dta[dta$X3=="exon",]$exon.id)
    exon.ranges.object <- GenomicRanges::reduce(exon.ranges.object)
    exon.ranges.object$exon=as.character(sprintf("%03d",1:length(exon.ranges.object)))
    ovl <- as_tibble(findOverlaps(exon.ranges.object,gene.ranges.object,type="any"))
    split.hits <- split(ovl$queryHits,ovl$subjectHits)
    exon.list <- lapply(names(split.hits),generate.exon.entry,a1=gene.ranges.object,a2=exon.ranges.object,split.hits=split.hits)
    exon.range <- do.call("c",exon.list)
    exon.range$exon=paste(exon.range$gene,exon.range$exon,sep="_")
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

## filter(gff,gene.id==filter(gff,exon.id=="ENSMUSE00000444479")$gene.id) %>%
##     slice.gene(.,window.shift=100,window.width=2e3) ->
##     aa

## aa <- slice.gene(filter(gff,gene.id=="ENSMUSG00000000001"),window.shift=100,window.width=2e3)

aa <- slice.gene(filter(gff,gene.id=="ENSMUSG00000033845"),window.shift=100,window.width=2e3)
aa <- slice.gene(filter(gff,gene.id=="ENSMUSG00000025933"),window.shift=100,window.width=2e3)
ll <- lineprof(slice.gene(filter(gff,gene.id=="ENSMUSG00000025933"),window.shift=100,window.width=2e3))

aa.old <- slice.gene(filter(gff,gene.id=="ENSMUSG00000025933"),window.shift=100,window.width=2e3)
aa.new <- slice.gene.new(filter(gff,gene.id=="ENSMUSG00000025933"),window.shift=100,window.width=2e3)

aa.old <- slice.gene(filter(gff,gene.id=="ENSMUSG00000033845"),window.shift=100,window.width=2e3)
aa.new <- slice.gene.new(filter(gff,gene.id=="ENSMUSG00000033845"),window.shift=100,window.width=2e3)



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
#    mutate(parent=paste("Parent=",parent,sep="")) %>%
#    unite("XX9",gene,exon,parent,sep=";") ->
    new.annotation.full
#    dplyr::select(X1,X2,

new.annotation.full[,] %>%
#    filter(X1%in% c("chr18","chr19")) %>%
    group_by(gene,exon,X1,X2,X3,X4,X5,X6,X7,X8) %>%
    do(X9=assemble.X9(.)) %>%
    tidy(X9) ->
    new.annotation.full.small

colnames(new.annotation.full.small)[ncol(new.annotation.full.small)] <- "X9"

new.annotation.out <- new.annotation.full.small[,paste("X",as.character(1:9),sep="")]
write_tsv(new.annotation.out,path="gencode_derived_sliding_window_genes.gff",col_names=F)


new.annotation.out %>%
    group_by(X1) %>%
    do(ll=(.)) %$%
    setNames(ll,X1) ->
    annotation.list

for (i in 1:length(annotation.list)){
    new.annotation.out <- annotation.list[[i]]
    write_tsv(new.annotation.out,path=paste("sliding_intron_",names(annotation.list)[i],".gff",sep=""),col_names=F)
}


aa <- slice.gene(filter(gff,gene.id=="ENSMUSG00000026127"))

new.annotation.full %>%
    separate(gene,into=c("original.gene","foo"),sep="_") %>%
    ungroup %>%
    mutate(foo=as.numeric(foo)) %>%
    mutate(ix=foo %% (2e3/100+1)) ->
    indexed.annotation

spacing <- 2e3/100

indexed.annotation %>%
    group_by(ix) %>%
    do(ll=(.)) %$%
    setNames(ll,ix) ->
    annotation.list


as.gr.spec <- function(x){
    GRanges(seqnames=x$X1,IRanges(start=x$X4,end=x$X5),strand=x$X7)
}

oo <- findOverlaps(as.gr.spec(annotation.list[[1]][annotation.list[[1]]$X3=="gene",]),as.gr.spec(annotation.list[[1]][annotation.list[[1]]$X3=="gene",]))
table(as_tibble(oo)[,1]-as_tibble(oo)[,2])

for (i in 1:length(annotation.list)){
    new.annotation.out <- annotation.list[[i]][,paste("X",as.character(1:9),sep="")]
    write_tsv(new.annotation.out,path=paste("gencode_derived_sliding_window_genes_",as.character(i),".gff",sep=""),col_names=F)
}

a


## dta <- filter(gff,gene.id==gff$gene.id[15])
## dta[dta$X3=="exon",]$X5[1] <- 7130615
## #if strand == "-" swap X4,X5
## gene.start <- dta[dta$X3=="gene",]$X4
## exons.start <- dta[dta$X3=="exon",]$X4
## exon.end <- max(dta[dta$X3=="gene",]$X5)
## gene.end <- exon.end
## dta[dta$X3=="gene",]$X5 <- gene.end
## borders <- seq(gene.start,gene.end,1e4)
## borders[length(borders)] <- gene.end
## sapply(exons.start,function(x){max(which(x>borders))})
## a1 <- IRanges(start=borders[1:(length(borders)-1)]+c(0,rep(1,(length(borders)-2))),end=borders[2:length(borders)])
## a2 <- IRanges(start=dta[dta$X3=="exon",]$X4,end=dta[dta$X3=="exon",]$X5)
## ovl <- as_tibble(findOverlaps(a2,a1,type="any"))

## split.hits <- split(ovl$subjectHits,ovl$queryHits)
## new.borders <- Reduce("c",lapply(split.hits,new.border,borders=a1))
## group_by(ovl,queryHits) %>%
##     mutate(lt=n()) %>%
##     ungroup %>%
##     filter(lt>1) ->
##     to.excise

## a1 <- sort(c(a1[-to.excise$subjectHits,],new.borders))

## new.border <- function(x,borders){
##     rv <- IRanges(start=c(),end=c())
##     if (length(x)>1){
##         new.border.start <- start(borders)[min(x)]
##         new.border.end <- end(borders)[max(x)]
##         rv <- IRanges(start=new.border.start,end=new.border.end)
##     }
##     rv
## }


## reset.borders <- function(x,borders){
##     if (length(x)>1){
##         new.border.start <- start(borders)[min(x)]
##         new.border.end <- end(borders)[max(x)]
##         borders <- c(borders[-x,],IRanges(start=new.border.start,end=new.border.end))
##     }
##     borders
## }


## dta <- filter(gff,gene.id==gff$gene.id[15])
## dta[dta$X3=="exon",]$X5[1] <- 7130615
## #if strand == "-" swap X4,X5
## gene.start <- dta[dta$X3=="gene",]$X4
## exons.start <- dta[dta$X3=="exon",]$X4
## exon.end <- max(dta[dta$X3=="gene",]$X5)
## gene.end <- exon.end
## dta[dta$X3=="gene",]$X5 <- gene.end
## borders <- seq(gene.start,gene.end,1e4)
## borders[length(borders)] <- gene.end
## sapply(exons.start,function(x){max(which(x>borders))})
## gene.names=paste(dta$gene.id[1],as.character(sprintf("%03d",1:(length(borders)-1))),sep="_")
## a1 <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=borders[1:(length(borders)-1)]+c(0,rep(1,(length(borders)-2))),end=borders[2:length(borders)]),gene=gene.names)
## a2 <- GRanges(seqnames="foo",strand=dta$X7[1],IRanges(start=dta[dta$X3=="exon",]$X4,end=dta[dta$X3=="exon",]$X5),exon=dta[dta$X3=="exon",]$exon.id)
## ovl <- as_tibble(findOverlaps(a2,a1,type="any"))
## split.hits <- split(ovl$queryHits,ovl$subjectHits)
## exon.list <- lapply(names(split.hits),function(x){rv <- intersect(a1[as.numeric(x),],a2[split.hits[[x]]]);rv$gene=a1[as.numeric(x),]$gene;rv$exon=a2[split.hits[[x]]]$exon;rv})
## exon.range <- do.call("c",exon.list)
## exon.range$exon=paste(exon.range$gene,exon.range$exon,sep="_")
## return.value <- as_tibble(c(a1,exon.range))[,c("start","end","gene","exon")]


bb <- sort(ceiling(100*runif(10)))

width=10
shift=5
cbind(bb[1]+shift*0:9,bb[1]+width+shift*0:9)
