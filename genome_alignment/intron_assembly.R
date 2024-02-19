library(GenomicRanges)
#library(seqQTL)
library(tidyverse)
#library(magrittr)

extract.ids <- function(string){
    gene.id <- str_extract(string,"ENSMUSG[0-9]*")
    exon.id <- str_extract(string,"ENSMUSE[0-9]*")
    transcript.id <- str_extract(string,"ENSMUST[0-9]*")
    list(gene.id=gene.id,exon.id=exon.id,transcript.id=transcript.id)
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
        id.string=paste("ID=",x$exon,";","exon_id=",true.gene.id,"_",x$exon,";","gene_id=",x$gene,";","Parent=",x$parent,sep="")
    }
                                        #    tibble(X9=id.string)
    id.string
}

rdc.helper <- function(to.reduce,reducer,ignore.strand=F){
    hits <- findOverlaps(to.reduce, reducer,ignore.strand=ignore.strand)
    grl <- extractList(reducer, as(hits, "List"))
    psd <- psetdiff(to.reduce, grl,ignore.strand=ignore.strand)
}



args = commandArgs(trailingOnly=TRUE)
gff.file=args[1]
no.cores=as.numeric(args[2])
#out.path=args[3]
#count.file="/extra/meisig/temp/slam_seq/results_exon_new/count_data/count_all_samples.csv"
#gff.file="/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3"
# no.cores=4

gff <- read_tsv(gff.file,col_names=F,guess_max=1e4,quote="",skip=7)

chr.names <- paste("chr",c(as.character(1:19),"X","Y","M"),sep="")

gff[,] %>%
    filter(X3 %in% c("gene","exon")) ->
    gff

ids <- mclapply(gff$X9,extract.ids,mc.cores=no.cores)

ids <- bind_rows(ids)
gff <- bind_cols(gff,ids)



gff %>%
    filter(X3=="gene") ->
    gff.gene

gff %>%
    filter(X3=="exon") ->
    gff.exon


gr.gff.gene <- GRanges(seqnames=gff.gene$X1,ranges=IRanges(start=gff.gene$X4,end=gff.gene$X5),strand=gff.gene$X7,gene=gff.gene$gene.id)
gr.gff.exon <- GRanges(seqnames=gff.exon$X1,ranges=IRanges(start=gff.exon$X4,end=gff.exon$X5),strand=gff.exon$X7)

hits <- findOverlaps(gr.gff.gene, gr.gff.exon,ignore.strand=T)
grl <- extractList(gr.gff.exon, as(hits, "List"))
#psetdiff(gr.gff.gene, grl,ignore.strand=T)

psd <- psetdiff(gr.gff.gene, grl,ignore.strand=T)
new.f <- mclapply(1:length(psd),function(x){y=psd[[x]];if (length(y)>0){y$gene=gr.gff.gene[x,]$gene};y},mc.cores=6)
new.f <- unlist(as(new.f[], "GRangesList"))

out.f <- new.f
out.f$gene <- make.unique(out.f$gene)

out.csv <- tibble(X1=as.character(seqnames(out.f)),X2="Havana",X3="gene",X4=start(out.f),X5=end(out.f),X6=".",X7=as.character(strand(out.f)),X8=".",X9=".",gene=out.f$gene)
transcripts <- out.csv
transcripts$parent=transcripts$gene
transcripts$X3="transcript"
exons <- out.csv
exons$parent=exons$gene
exons$X3="exon"
exons$exon=paste("E",1:nrow(exons),sep="_")
exons$parent <- paste(exons[,]$parent,"_T",sep="")

out.csv <- bind_rows(out.csv,transcripts,exons)
#out.csv <- bind_rows(out.csv)
out.csv <- out.csv[order(out.csv$gene),]

assembled.X9 <- unlist(mclapply(1:nrow(out.csv),function(x){assemble.X9(out.csv[x,])},mc.cores=6))
out.csv$X9=assembled.X9

write_tsv(out.csv[,1:9],path="introns_only.gff",col_names=F)

## out.f$gene=str_replace(out.f$gene[],"\\.[0-9]+","")
    

## out.f %>%
##     as_tibble() %>%
##     group_by(gene) %>%
##     group_split() ->
##     intron.gff

## out.f %>%
##     as_tibble() %>%
##     group_by(gene) %>%
##     group_keys() ->
##     intron.gff.annot

## names(intron.gff) <- intron.gff.annot$gene


## save(intron.gff,file="intron_mask.RData")


##     gene.ranges.object <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=as.numeric(borders[,1]),end=as.numeric(borders[,2])),gene=gene.names,type="gene")
##     transcripts <- GRanges(seqnames="foo",strand=dta$X7[1],ranges=IRanges(start=as.numeric(borders[,1]),end=as.numeric(borders[,2])),gene=gene.names,parent=gene.names,type="transcript")

##     exon.ranges.object <- GRanges(seqnames="foo",strand=dta$X7[1],IRanges(start=dta[dta$X3=="exon",]$X4,end=dta[dta$X3=="exon",]$X5),exon=dta[dta$X3=="exon",]$exon.id)


## library(GenomicFeatures)

## gff_file <- system.file("extdata", "GFF3_files", "/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3",package="GenomicFeatures")
## txdb <- makeTxDbFromGFF("/extra/meisig/temp/slam_seq/reference/gencode.vM14.basic.annotation.gff3", format="gff3")
## txdb


## trak2_txs <- transcriptsBy(txdb, by="gene")
## trak2_txs

## trak2_inbytx <- intronsByTranscript(txdb, use.names=TRUE)

## lapply(1:length(exons),function(x){intersect(gene,exons[x,])})
