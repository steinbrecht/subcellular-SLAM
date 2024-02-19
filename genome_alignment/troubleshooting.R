new.annotation %>%
    group_by(gene) %>%
    summarise(st=min(start),ed=max(end)) ->
    start.end

new.annotation %>%
    filter(type=="transcript") %>%
    group_by(gene) %>%
    summarise(st=min(start),ed=max(end)) ->
    start.end.trans

max(start.end$st-start.end.trans$st)

min(start.end$ed-start.end.trans$ed)

chr5 <- bind_cols(annotation.list[["chr5"]],bind_rows(lapply(annotation.list[["chr5"]]$X9,extract.ids)))

chr5 %>%
    group_by(gene.id) %>%
    group_split() ->
    chr5.split

lapply(chr5.split,nrow)


for (i in 1:length(chr5.split)){
    new.annotation.out <- chr5.split[[i]][,paste("X",as.character(1:9),sep="")]
    write_tsv(new.annotation.out,path=paste(out.path,"/troubleshoot_",chr5.split[[i]]$gene.id[1],".gff",sep=""),col_names=F)
}




chr5 %>% filter(gene.id=="ENSMUSG00000034842") ->foo

foo %>% filter(grepl("_gene",X9)) -> bar

autoplot(GRanges(seqnames=bar$X1,IRanges(start=bar$X4,end=bar$X5)))

gff %>% filter(exon.id=="ENSMUSE00000542411")  %>% as.data.frame
