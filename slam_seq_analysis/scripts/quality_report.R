args = commandArgs(trailingOnly=TRUE)
data.dir=args[1]
mut.table.dir=args[2]
align.data.dir=args[3]
script.dir=args[4]

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

conv.rates <- read_csv(paste(data.dir,"/conversion_probability_fit_mixture_model.csv",sep="")) %>%
    filter(!grepl("IAA",sample)) %>%
    filter(!grepl("ameres",fraction)) %>%
    convert.sample %>%
    mutate(pr=p.conv)


conv.rates %>%
    filter(sample!="pc") %>%
    filter(time>0) %>%
    mutate(rep=as.character(rep)) %>%
    filter(rep%in% as.character(1:4)) %>%
    ggplot(aes(x=time,y=pr,group=interaction(sample,fraction,rep),colour=fraction))+
    geom_line()+
    geom_point()+ 
    facet_wrap(~rep)+
    scale_colour_brewer(palette="Set1")+
    scale_x_continuous("Time [h]")+
    scale_y_continuous("Conversion probability")->
    p1

save_plot(filename=paste(data.dir,"/conversion_probability_different_fractions.pdf",sep=""),plot=p1,base_height=5,base_aspect_ratio=2)


conv.rates %>%
    filter(fraction=="chase") %>%
    mutate(rep=as.character(rep)) %>%
    ggplot(aes(x=time,y=pr,group=interaction(sample,fraction,rep),colour=sample,linetype=rep))+
    geom_line()+
    geom_point()+ 
    scale_colour_brewer(palette="Set1")+
    scale_x_continuous("Time [h]")+
    scale_y_continuous("Conversion probability")->
    p1

save_plot(filename=paste(data.dir,"/conversion_probability_only_chase.pdf",sep=""),plot=p1,base_height=5,base_aspect_ratio=2)

conv.rates %>%
    filter(fraction=="chase") %>%
    filter(time>0) %>%
    mutate(rep=as.character(rep)) %>%
#    filter(sample %in% c("n","t","m","c")) %>%
    filter(rep%in% as.character(1:4)) %>%
    mutate(log.fraction.labeled=log(a)) %>%
    ggplot(aes(x=time,y=log.fraction.labeled,group=interaction(sample,fraction,rep),colour=sample,linetype=rep))+
    geom_line()+
    geom_point()+ 
    scale_x_continuous("Time [min]")->
    p1

save_plot(filename=paste(data.dir,"/mixture_parameter_only_chase.pdf",sep=""),plot=p1,base_height=5,base_aspect_ratio=2)


conv.rates %>%
    filter(fraction!="chase") %>%
    filter(sample!="pc") %>%
    filter(time>0) %>%
    filter(rep%in% as.character(1:4)) %>%
    mutate(fraction.labeled=(a)) %>%
    ggplot(aes(x=time,y=fraction.labeled,group=interaction(sample,fraction,rep),colour=fraction,linetype=rep))+
    geom_line()+
    geom_point()+ 
    scale_x_continuous("Time [min]")+
    scale_colour_brewer(palette="Set1")->
    p1

save_plot(filename=paste(data.dir,"/mixture_parameter_different_fractions.pdf",sep=""),plot=p1,base_height=5,base_aspect_ratio=2)


read.multiple.files <- function(path,read.function,pattern,col.types=NULL,no.cores=1){
    all.files <- list.files(path=path,pattern=pattern,full.names=TRUE)
    names.all.files <- list.files(path=path,pattern=pattern,full.names=FALSE)
    all.data <- mclapply(all.files,read.function,col_names=FALSE,col_types=col.types,mc.cores=no.cores)
    all.data <- bind_rows(all.data,.id="sample")
    all.data %>%
        mutate(sample=names.all.files[as.numeric(sample)]) ->
        all.data
}

frequency.data <- read.multiple.files(path=mut.table.dir,read_tsv,pattern="freq.*tsv",col.types=c(.default = col_integer()),no.cores=6)

frequency.data %>%
    mutate(sample=str_replace(sample,"frequencies_","")) %>%
    mutate(sample=str_replace(sample,".tsv","")) %>%
    separate(sample,c("type","sample","fraction"),sep="[.]") ->
    frequency.data


frequency.data %>%
    filter(X1>0) %>%
    filter(type=="A2G") %>%
    select(sample,fraction,X1,X2,X3) %>%
    group_by(sample,fraction,X1) %>%
    mutate(p.error=ml.estimate.pconv(data=data_frame(N=X1,k=X2,akn=X3))) %>%
    select(sample,fraction,X1,p.error) %>%
    distinct ->
    error.estimates


error.estimates %>%
    filter(sample=="n0min1") %>%
    filter(X1<100) %>%
    ggplot(aes(x=X1,y=p.error))+
    geom_line()+
    scale_y_log10() ->
    p1

save_plot(filename=paste(data.dir,"/dependence_error_rate_no_ts.pdf",sep=""),plot=p1,base_height=4,base_aspect_ratio=2)


read_csv(paste(data.dir,"/conversion_probability_fit_mixture_model.csv",sep="")) %>%
    mutate(b=1-a) %>%
    select(sample,fraction,a,b,cost,p.conv) ->
    conv.rates

    

frequency.data %>%
    left_join(error.estimates,by=c("sample","fraction","X1")) %>%
    left_join(conv.rates,by=c("sample","fraction")) %>%
    select(type,sample,fraction,N=X1,k=X2,akn=X4,a,b,cost,p.conv,p.error) %>%
    filter(type=="T2C") %>%
    group_by(sample,N) %>%
    mutate(sum.akn=sum(akn)) %>%
    filter(N %in% c(35,seq(10,60,5))) %>%
#    filter(N==35) %>%
    group_by(sample,fraction,N,k,akn) %>%
    do(data_frame(fit=sol.binomial.mixt(.$a[1],.$b[1],.$p.conv[1],.$p.error[1],.))) %>%
    unnest %>%
    gather("type","value",akn:fit) ->
    fit.mixture.out


# fit.mixture.out %>%
#     ungroup %>%
#     filter(fraction=="chase") %>%
#     distinct %>%
#     convert.sample %>%
#     filter(time==360,sample=="n") %>%
#     filter(k<12) %>%
#     mutate(k=as.factor(k)) %>%
#     mutate(events=value) %>%
#     mutate(events=log10(events+1)) %>%
#     mutate(type=ifelse(type=="akn","data","fit")) %>%
#     ggplot(aes(x=k,y=events,colour=type,group=type))+
#     geom_line()+
#     facet_wrap(time~rep+sample+N,scales="free_y")+
#     scale_y_continuous("Log10 # events +1")+
#     scale_x_discrete("k")->
#     p1
            

# save_plot(filename=paste(data.dir,"/fit_pconv_chase_360.pdf",sep=""),plot=p1,base_height=7,base_aspect_ratio=2)

write_csv(fit.mixture.out,path=paste(data.dir,"/frequency_fit_binomial_mixture_model.csv",sep=""))

fit.mixture.out %>%
    ungroup %>%
    filter(fraction=="nuc") %>%
    distinct %>%
    convert.sample %>%
    filter(time==180,sample=="n") %>%
    filter(k<12) %>%
    mutate(k=as.factor(k)) %>%
    mutate(events=value) %>%
    mutate(events=log10(events+1)) %>%
    mutate(type=ifelse(type=="akn","data","fit")) %>%
    ggplot(aes(x=k,y=events,colour=type,group=type))+
    geom_line()+
    facet_wrap(time~rep+sample+N,scales="free_y")+
    scale_y_continuous("Log10 # events +1")+
    scale_x_discrete("k")->
    p1
            

save_plot(filename=paste(data.dir,"/fit_pconv_nuc_180.pdf",sep=""),plot=p1,base_height=7,base_aspect_ratio=2)

fit.mixture.out %>%
    ungroup %>%
    filter(fraction=="cyto") %>%
    distinct %>%
    convert.sample %>%
    filter(time==180) %>%
    filter(k<12) %>%
    mutate(k=as.factor(k)) %>%
    mutate(events=value) %>%
    mutate(events=log10(events+1)) %>%
    mutate(type=ifelse(type=="akn","data","fit")) %>%
    ggplot(aes(x=k,y=events,colour=type,group=type))+
    geom_line()+
    facet_wrap(time~rep+sample+N,scales="free_y")+
    scale_y_continuous("Log10 # events +1")+
    scale_x_discrete("k")->
    p1
            

save_plot(filename=paste(data.dir,"/fit_pconv_cyto_180.pdf",sep=""),plot=p1,base_height=7,base_aspect_ratio=2)



gene.data=read_csv(paste(data.dir,"/gene_t2c_over_t_normalized.csv",sep=""))

gene.data %>%
    ungroup %>%
    mutate(rep=as.character(rep)) %>%
    filter(gene_biotype=="protein_coding") %>%
#    filter(T.intron>5e2) %>%
    filter(time!=0) %>%
    gather("type","t2c.rate",t2c.over.t.intron,t2c.over.t.exon) %>%
    mutate(type=ifelse(type=="t2c.over.t.intron","intron","exon")) %>%
    unite("fraction",sample,fraction,type) %>%
    filter(fraction %in% c("t_total_intron","n_nuc_intron","c_cyto_intron","m_mem_intron","pH_poly_intron","pL_poly_intron","n_chase_intron","m_chase_intron","c_chase_intron")) %>%
    select(symbol,gene,time,rep,fraction,t2c.rate) %>%
    distinct %>%
    ungroup %>%
    ggplot(aes(x=t2c.rate,group=interaction(time,fraction,rep),linetype=rep,colour=time))+
    geom_vline(xintercept=1,colour="darkgray")+
    stat_ecdf()+
    coord_cartesian(xlim=c(0,2))+
    scale_colour_gradient(low="darkblue",high="gold")+
    facet_wrap(~fraction,scales="free_y")+
    scale_x_continuous("Fraction new mRNA")->
    p1

save_plot(filename=paste(data.dir,"/ecdf_replicates.pdf",sep=""),plot=p1,base_height=7,base_aspect_ratio=2)

gene.data %>%
    ungroup %>%
    mutate(rep=as.character(rep)) %>%
    filter(gene_biotype=="protein_coding") %>%
#    filter(T.intron>5e2) %>%
    filter(time!=0) %>%
    gather("type","t2c.rate",t2c.over.t.intron,t2c.over.t.exon) %>%
    mutate(type=ifelse(type=="t2c.over.t.intron","intron","exon")) %>%
    unite("fraction",sample,fraction,type) %>%
    filter(fraction %in% c("t_total_exon","n_nuc_exon","c_cyto_exon","m_mem_exon","pH_poly_exon","pL_poly_exon","n_chase_exon","m_chase_exon","c_chase_exon")) %>%
    select(symbol,gene,time,rep,fraction,t2c.rate) %>%
    distinct %>%
    ungroup %>%
    ggplot(aes(x=t2c.rate,group=interaction(time,fraction,rep),linetype=rep,colour=time))+
    geom_vline(xintercept=1,colour="darkgray")+
    stat_ecdf()+
    coord_cartesian(xlim=c(0,2))+
    scale_colour_gradient(low="darkblue",high="gold")+
    facet_wrap(~fraction,scales="free_y")+
    scale_x_continuous("Fraction new mRNA")->
    p1

save_plot(filename=paste(data.dir,"/ecdf_replicates_exon.pdf",sep=""),plot=p1,base_height=7,base_aspect_ratio=2)


half.lifes.all <- read_csv(paste(script.dir,"/half_lifes_literature.csv",sep=""))
colnames(half.lifes.all)[colnames(half.lifes.all)=="ensembl_gene_id"]="gene"

filter(half.lifes.all,gene_biotype=="protein_coding",variable=="herzog_half.life") %>%
    filter(value>6) %$%
    unique(symbol) ->
    long.lived.genes

gene.data %>%
    filter(gene_biotype=="protein_coding") %>%
    group_by(gene) %>%
    mutate(mn=mean(T.exon,na.rm=TRUE)) %>%
    ungroup %>%
    filter(time>0) %>%
    mutate(mn=ntile(mn,10)) %>%
    filter(mn>=10) %>%
 #   filter(fraction%in%c("nuc","poly")) %>%
#    filter(fraction=="nuc"|sample=="pL"|sample=="pH") %>%
    filter(symbol%in% long.lived.genes) %>%
    filter(gene %in% sample(unique(gene),30)) %>%
    ggplot(aes(x=time,y=t2c.over.t.exon,group=interaction(rep,sample,symbol),colour=interaction(fraction,sample),linetype=as.factor(rep)))+
    geom_line()+
    geom_point()+
    facet_wrap(~symbol,scale="free_y")->
    p1

save_plot(filename=paste(data.dir,"/kinetics_long_lived_examples.pdf",sep=""),plot=p1,base_height=7,base_aspect_ratio=2)




# gene.data %>%
#     filter(gene_biotype=="protein_coding") %>%
#     filter(fraction=="chase") %>%
#     filter(t2c.over.t.exon>0) %>%
# #    filter(time>0 & time <300) %>%
#     filter(time>0) %>%
#     mutate(t2c.over.t.exon=log(t2c.over.t.exon)) %>%
#     group_by(sample,gene,rep) %>%
#     do(as_data_frame(t(as.matrix(lm(t2c.over.t.exon~time,data=.,na.action=NULL)$coefficients))))  ->
#     fits.chase


# write_csv(fits.chase,path=paste(data.dir,"/fits_chase_data.csv",sep=""))

# fits.chase %>%
#     ungroup %>%
#     left_join(half.lifes.all,by="gene")%>%
#     filter(variable=="herzog_half.life") %>%
#     mutate(rep=as.character(rep)) %>%
#     group_by(variable,sample,rep) %>%
#     mutate(thalf=-log(2)/time) %>%
#     filter(thalf>0) %>%
#     filter(value>0) %>%
#     summarise(cor=cor(value,thalf,method="spearman")) %>%
#     ungroup %>%
#     arrange(cor) ->
#     correlations.chase

# write_csv(correlations.chase,path=paste(data.dir,"/correlations_literature_half_lifes.csv",sep=""))


# fits.chase %>%
#     ungroup %>%
# #    filter(sample=="m") %>%
#     left_join(half.lifes.all,by="gene")%>%
#     filter(variable=="herzog_half.life") %>%
#     mutate(rep=as.character(rep)) %>%
#     mutate(value=value*60) %>%
#     group_by(variable,sample,rep) %>%
#     mutate(thalf.est=-log(2)/time) %>%
#     filter(thalf.est>0) %>%
#     filter(value>0) %>%
#     ggplot(aes(x=log10(value),y=log10(thalf.est)))+
#     geom_point(alpha=0.1,size=0.1)+
#     geom_abline(intercept=0,slope=1,colour="red",alpha=0.5)+
#     facet_wrap(sample~rep,ncol=4,scales="free_y")+
#     scale_x_continuous("Log10 literature half life [min]")+
#     scale_y_continuous("Log10 estimated half life [min]")->
#     p1

# save_plot(filename=paste(data.dir,"/correlation_chase_literature_new.pdf",sep=""),plot=p1,base_height=6,base_aspect_ratio=1.5)


read.multiple.files <- function(path,read.function,pattern,col.types=NULL,no.cores=1){
    all.files <- list.files(path=path,pattern=pattern,full.names=TRUE)
    names.all.files <- list.files(path=path,pattern=pattern,full.names=FALSE)
    all.data <- mclapply(all.files,read_tsv,col_names=FALSE,col_types=col.types,mc.cores=no.cores)
    all.data <- bind_rows(all.data,.id="sample")
    all.data %>%
        mutate(sample=names.all.files[as.numeric(sample)]) ->
        all.data
}


dd <- read.multiple.files(align.data.dir,read.function="read_tsv","*.Log.final.out")

dd %>%
                                        #    filter(grepl("mismatch",X1)|grepl("mismatch",X1)) %>%
#    filter(grepl("\\%",X1)|grepl("Average",X1)) %>%
    filter(grepl("\\%",X1)) %>%
    separate(sample,into=c("sample","foo","read"),sep="\\.") %>%
    separate(sample,into=c("foo1","sample","foo2","foo3"),sep="_") %>%
    select(sample,X1,X2) %>%
    convert.sample ->
    dd


dd %>%
    filter(sample %in% c("t","n","m","c","pL","pH")) %>%
    mutate(nn=as.numeric(str_extract(X2,"[0-9]*\\.*[0-9]*"))) %>%
    filter(!grepl("chimeric",X1)) %>%
    ggplot(aes(x=time,y=nn,group=interaction(sample,rep),colour=rep))+
    geom_line()+
    geom_point()+
    facet_grid(X1~sample,scales="free_y") ->
    p1


save_plot(file=paste(data.dir,"/star_log.pdf",sep=""),plot=p1,base_aspect_ratio=0.5,base_height=15)

