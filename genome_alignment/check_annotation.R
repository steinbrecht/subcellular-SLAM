
extract.ids <- function(string){
    gene.id <- str_extract(string,"ENSMUSG[0-9]+_[0-9]+|ENSMUSG[0-9]+_ENSMUSE[0-9]+|ENSMUSG[0-9]+_gene")
    exon.id <- str_extract(string,"ENSMUSE[0-9]*")
    transcript.id <- str_extract(string,"ENSMUST[0-9]*")
    list(gene.id=gene.id,exon.id=exon.id,transcript.id=transcript.id)
}


chr.names <- paste("chr",c(as.character(1:19),"X","Y","M"),sep="")

chr <- list()
for (i in chr.names){
    print(i)
    chr[[i]] <- read_tsv(paste("/extra/meisig/temp/slam_seq/reference/RSEM_SNAKE/sliding_intron_parallel_",i,".gff",sep=""),col_names=F)
}

id.addition <- function(annot.frame){
    ids <- mclapply(annot.frame$X9,extract.ids,mc.cores=6)
    ids <- bind_rows(ids)
    annot.frame <- bind_cols(annot.frame,ids)
}


chr <- lapply(chr,id.addition)
chr <- bind_rows(chr)
chr %>%
    separate(gene.id,into=c("gene","sub.id"),sep="_") %>%
    mutate(sub.id.no=as.numeric(sub.id)) ->
    chr


chr %>%
    select(gene,X4,X5,X7,gene,sub.id,sub.id.no,exon.id) ->
    chr

chr %>%
    filter(sub.id.no==min(sub.id.no)|sub.id.no==max(sub.id.no)|is.na(sub.id.no)) %>%
    group_by(gene) %>%
    mutate(end.wn=ifelse(X7=="+",max(X5[!is.na(sub.id.no)]),min(X4[!is.na(sub.id.no)]))) %>%
    mutate(end.ex=ifelse(X7=="+",max(X5[is.na(sub.id.no)]),min(X4[is.na(sub.id.no)]))) ->
    foo


sum(abs(foo$end.wn-foo$end.ex)==0)


chr %>%
    filter(X1=="chr2") %>%
#    filter(!is.na(exon.id)) %>%
    summarise(sum(duplicated(sub.id)),sum(duplicated(exon.id)))

chr %>%
    filter(X1=="chr2") %>%
    group_by(gene) %>%
    mutate(mn=min(X4[X3=="gene"]),mx=max(X5[X3=="gene"])) %>%
    select(X3,X4,X5,mn,mx,gene,sub.id,exon.id) %>%
    mutate(X4=X4-mn,X5=mx-X5) %>%
    filter(X5<0)


tr <- read_tsv(paste("/extra/meisig/temp/slam_seq/reference/RSEM_SNAKE/trouble3.gff"),col_names=F)

View(tr[duplicated(paste(tr$X3,tr$X4,tr$X5,sep="_")),])

tr %>%
    filter(X3=="exon") %>%
    select(-X1,-X2,-X6,-X8) %>%
    as.data.frame

tr %>%
    filter(X3=="exon") %>%
    select(-X1,-X2,-X6,-X8) %>%
    filter(X4==32774391) %>%
    as.data.frame
