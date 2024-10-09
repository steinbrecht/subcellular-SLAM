library(dplyr)
library(readr)
library(tidyr)
library(lhs)
                                        #library(purrr)
library(magrittr)
#library(stringr)
#library(superslam)
library(parallel)


args = commandArgs(trailingOnly=TRUE)
in.file.t2c=args[1]
in.file.a2g=args[2]
out.file=args[3]

no.random=as.numeric(as.character(args[4]))
random.fact=as.numeric(as.character(args[5]))
no.cores=as.numeric(as.character(args[6]))
use.intronic.reads=args[7]

source("scripts/mixture_optimization.R")


filter.fun.intronic=function(x){
  x %>%    
    select(num_Ts,num_T2Cs,freq_introns) %>%
    filter(num_Ts>0) %>%
    rename(freq_exons=freq_introns)
}

filter.fun.exonic=function(x){
  x %>%    
    select(num_Ts,num_T2Cs,freq_exons) %>%
    filter(num_Ts>0)
}


if (use.intronic.reads){
   filter.fun=filter.fun.intronic
}else{
   filter.fun=filter.fun.exonic
}


frequency.t2c <- read_tsv(in.file.t2c,col_types=c(.default = col_integer()),col_names=F)
frequency.a2g <- read_tsv(in.file.a2g,col_types=c(.default = col_integer()),col_names=F)
cn=c("num_Ts","num_T2Cs","freq_exons","freq_introns")
colnames(frequency.t2c)=cn
colnames(frequency.a2g)=cn



frequency.a2g %>%
    select(num_Ts,num_T2Cs,freq_exons) %>%
    rename(N=num_Ts,k=num_T2Cs,akn=freq_exons) ->
    all.errors


frequency.t2c %>%
    filter.fun ->
    all.events

mixture.fit <- multistep.per.T.error(data=all.events,data.error=all.errors,no.random=no.random,threshold.excluded.events=25,random.fact=random.fact)

write_csv(mixture.fit,path=out.file)
