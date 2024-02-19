library(stringr)
library(parallel)
library(tidyverse)



contamination.matrix=function(f.vec){
    fmc=f.vec[1]
    fnc=f.vec[2]
    fcm=f.vec[3]
    fnm=f.vec[4]
    fcn=f.vec[5]
    fmn=f.vec[6]
    row.1=c(1,fmc,fnc)/(1+fmc+fnc)
    row.2=c(fcm,1,fnm)/(1+fcm+fnm)
    row.3=c(fcn,fmn,1)/(1+fcn+fmn)
    rbind(row.1,row.2,row.3)
}

cost.function=function(x,fit.data,true.counts){
    prediction=contamination.matrix(x) %*% true.counts
    prediction=log10(prediction)
    fit.data=log10(fit.data)
    sum((fit.data-prediction)^2)
}


opt.wrapper <- function(data,rep.reference,rep.fit){
#    reference.data=data[[rep.reference]]
    ref.hits=lapply(rep.reference,function(x){grepl(x,names(fit.data.in))})
    reference.data=lapply(ref.hits,function(x){Reduce("+",fit.data.in[x])/sum(x)})
    reference.data=rbind(reference.data[[1]][1,],reference.data[[2]][2,],reference.data[[3]][3,])
    fit.data=data[[rep.fit]]
    f.start=rep(1e-3,6)
    optim(par=f.start,fn=cost.function,upper=rep(1,length(f.start)),lower=rep(0,length(f.start)),method="L-BFGS-B",fit.data=fit.data,true.counts=reference.data)
}

args = commandArgs(trailingOnly=TRUE)
fit.data=args[1]
no.cores=as.numeric(as.character(args[2]))
output.file=args[3]

load(fit.data)

reference.source=1:6

expand.grid(ref.1=paste(as.character(reference.source),"_",sep=""),ref.2=paste(as.character(reference.source),"_",sep=""),ref.3=paste(as.character(reference.source),"_",sep=""),fit=names(fit.data.in),stringsAsFactors=F) %>%
    as_tibble ->
    todo


opt.out=mclapply(1:nrow(todo),function(x){print(x/nrow(todo));opt.wrapper(data=fit.data.in,rep.reference=as.character(todo[x,c("ref.1","ref.2","ref.3")]),rep.fit=as.character(todo[x,"fit"]))},mc.cores=no.cores)


bind_cols(todo,enframe(unlist(lapply(opt.out,"[[","value")))) %>%
    group_by(ref.1,ref.2,ref.3) %>%
    summarise(md=median(value)) %>%
    arrange(md) ->
    cost.frame

out.pars=do.call("rbind",lapply(opt.out,"[[","par"))
colnames(out.pars)=c("fmc","fnc","fcm","fnm","fcn","fmn")

out.pars %>%
    as_tibble %>%
    bind_cols(todo) ->
    out.pars

out.pars %>%
    left_join(cost.frame) %>%
    rename(cost=md) %>%
    arrange(cost) ->
    out.pars

write_csv(x=out.pars,output.file)
