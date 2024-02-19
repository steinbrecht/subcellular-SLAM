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

retain.exon.info <- function(transcript.info,dta){
    gene.start=transcript.info$X4
    gene.end=transcript.info$X5
    transcript.id=transcript.info$transcript.id
    gene=paste(dta$gene.id[1],"_gene",sep="")
    gene.tbl=tibble(type="gene",parent=NA,start=gene.start,end=gene.end,gene=gene,exon=NA)
    transcript.tbl=tibble(type="transcript",parent=gene,start=gene.start,end=gene.end,gene=gene,exon=NA)
#    exon.tbl <- tibble(type="exon",parent=paste(dta$gene.id[1],transcript.id,sep="_"),start=dta$start,end=dta$end,gene=gene,exon=dta$exon.id,count=NA)
    exon.tbl <- tibble(type="exon",parent=paste(dta$gene.id[1],"T",sep="_"),start=dta$start,end=dta$end,gene=gene,exon=dta$exon.id,count=NA)
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
    psd
}


slice.gene.new <- function(info.list,window.width,window.shift){
    dta=info.list[["original"]]
    intron.restriction=info.list[["intron.windows"]]
    gene.start <- min(dta[dta$X3=="exon",]$X4)
    gene.end <- max(dta[dta$X3=="exon",]$X5)
    if (gene.end-gene.start<window.width){
#        borders=cbind(gene.start,gene.end)
        borders=cbind(intron.restriction$start,intron.restriction$end)
    }else{
        no.shifts <- floor((gene.end-gene.start-window.width)/window.shift)
        borders <- cbind(gene.start+window.shift*0:no.shifts,gene.start+window.width+window.shift*0:no.shifts)
        borders[nrow(borders),2] <- gene.end
    }
    gene.names=paste(dta$gene.id[1],as.character(sprintf("%05d",1:nrow(borders))),sep="_")
    gene.ranges.object <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=as.numeric(borders[,1]),end=as.numeric(borders[,2])),gene=gene.names,type="gene")
    ab <- GRanges(seqnames="foo",ranges=IRanges(start=intron.restriction$start,end=intron.restriction$end))
    gene.ranges.object <- rdc.helper(gene.ranges.object,ab)
#    gene.ranges.object=gene.ranges.object[width(gene.ranges.object)>window.width/4,]
    if (length(gene.ranges.object)>0){
    gene.names=paste(dta$gene.id[1],as.character(sprintf("%05d",1:length(gene.ranges.object))),sep="_")
    gene.ranges.object$gene=gene.names
        transcripts <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=start(gene.ranges.object),end=end(gene.ranges.object)),gene=gene.ranges.object$gene,parent=gene.ranges.object$gene,type="transcript")
    }else{
        transcripts <- GRanges(seqnames=NULL,strand=NULL,ranges=IRanges(start=NULL,end=NULL),gene=NULL,parent=NULL,type=NULL)
        gene.ranges.object <- GRanges(seqnames=NULL,strand=NULL,ranges=IRanges(start=NULL,end=NULL),gene=NULL,parent=NULL,type=NULL)
    }
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

restrict.fn <- function(example.gff){
    expressed.transcript <- unique(example.gff[!is.na(example.gff$TPM),]$transcript.id)
    example.gff %>%
        filter(transcript.id==expressed.transcript) %>%
        filter(X3=="transcript")%>%
        select(X4,X5) %>%
        as.numeric ->
        borders
    example.gr <- GRanges(seqnames=example.gff$X1,ranges=IRanges(start=example.gff$X4,end=example.gff$X5),strand=example.gff$X7,gene.id=example.gff$gene.id,exon.id=example.gff$exon.id,transcript.id=example.gff$transcript.id,type=example.gff$X3)
    example.gr <- restrict(example.gr,start=borders[1],end=borders[2])
    exonic <- example.gr[!is.na(example.gr$exon.id)&example.gr$type=="exon",]
    intronic <- gaps(exonic)[-1,]
    exonic <- exonic[!duplicated(exonic$exon.id),]
    exonic <- unique(exonic)
    example.gff %>%
        filter(transcript.id==expressed.transcript) ->
        example.gff
    list(gff=example.gff,exonic=exonic,intronic=intronic)
}


extract.info <- function(gff.slice,intron.ranges){
    print(gff.slice$gene.id[1])
    out.d <- restrict.fn(gff.slice)
    if (length(out.d[[3]])>0){
        out.d[[3]] <- rdc.helper(to.reduce=out.d[[3]],reducer=intron.ranges,ignore.strand=F)
    }
    out.d[[3]] <- as_tibble(out.d[[3]])
    if (nrow(out.d[[3]])>0){
        aa=slice.gene.new(list(original=out.d[[1]],intron.windows=out.d[[3]]),window.width=4e3,window.shift=4e2)
    }else{
        aa=tibble(type=character(),parent=character(),start=double(),end=double(),gene=character())
    }
    transcript.info=out.d[[1]][out.d[[1]]$X3=="transcript",]
    ab=retain.exon.info(transcript.info,as_tibble(out.d[[2]]))
    bind_rows(aa,ab)
}


args = commandArgs(trailingOnly=TRUE)
print(args)
gff.file=args[1]
no.cores=as.numeric(as.character(args[2]))
out.path=args[3]
#gff.file="/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3"
#no.cores=4
print(no.cores)

#gff.intron=read_tsv("introns_only.gff",col_names=F,guess_max=1e4,quote="",skip=7) %>%
#    filter(X3=="gene")
#ids.intron <- mclapply(gff.intron$X9,extract.ids,mc.cores=no.cores)
#ids.intron <- bind_rows(ids.intron)
#gff.intron <- bind_cols(gff.intron,ids.intron)
#intron.ranges <- GRanges(seqnames=gff.intron$X1,strand=gff.intron$X7,ranges=IRanges(start=gff.intron$X4,end=gff.intron$X5))


gff <- read_tsv(gff.file,col_names=F,guess_max=1e4,quote="",skip=7)
gff.exon <- filter(gff,X3=="exon")
exon.ranges <- GRanges(seqnames=gff.exon$X1,strand=gff.exon$X7,ranges=IRanges(start=gff.exon$X4,end=gff.exon$X5))
intron.ranges <- gaps(exon.ranges)

chr.names <- paste("chr",c(as.character(1:19),"X","Y","M"),sep="")

ids <- mclapply(gff$X9,extract.ids,mc.cores=no.cores)
ids <- bind_rows(ids)
gff <- bind_cols(gff,ids)
gff <- gff[!is.na(gff$gene.id),]

rsem.top=read_csv("rsem_top_expressed_transcripts_total.csv") %>%
    rename(transcript.id=transcript,gene.id=gene)


gff %>%
    left_join(rsem.top,by=c("gene.id","transcript.id")) %>%
    filter(gene.id %in% gene.id[!is.na(TPM)]) ->
    gff


gff %>%
    group_by(gene.id) %>%
    group_split() ->
    split.gff

gff %>%
    group_by(gene.id) %>%
    group_keys() ->
    split.gff.annot


names(split.gff) <- split.gff.annot$gene.id


                                        #extract.info(split.gff[["ENSMUSG00000000295"]])
#extract.info(split.gff[["ENSMUSG00000019836"]],intron.ranges=intron.ranges)
#extract.info(split.gff[["ENSMUSG00000039021"]],intron.ranges=intron.ranges)

#new.annotation <- mclapply(split.gff[sample(names(split.gff),500)],extract.info,mc.cores=6)
new.annotation <- mclapply(split.gff[],extract.info,intron.ranges=intron.ranges,mc.cores=no.cores)
new.annotation <- bind_rows(new.annotation)

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
