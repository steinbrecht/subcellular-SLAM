source("scripts/superslam_algorithm.R")

convert.sample <- function(x){
    x%>%
    mutate(sample=str_replace(sample,"min","_")) %>%
    # mutate(sample=str_replace(sample,"ch_","")) %>%
    separate(sample,c("sample","rep"),sep="_") %>%
    mutate(time=str_extract(sample,"(\\d)+")) %>%
    mutate(sample=str_replace(sample,"[0-9].","")) %>%
    mutate(sample=str_replace(sample,"[0-9]","")) %>%
    mutate(time=as.numeric(time))
}

sol.binomial.mixt <- function(a,b,p.conv,p.error,data){
    pred.conv=data$sum.akn*dbinom(x=data$k,size=data$N,prob=p.conv)
    pred.error=data$sum.akn*dbinom(x=data$k,size=data$N,prob=p.error)
    pred=a/(a+b)*pred.conv+b/(a+b)*pred.error
}


cost.fn.binomial <- function(params,p.conv,p.error,data){
    a=params[1]
    b=1-a
    pred=sol.binomial.mixt(a,b,p.conv,p.error,data)
    res=data$akn
    sum(data$akn/(sum(data$akn))*(pred-res)^2)
}


binomial.solution.wrapper <- function(a,b,p.conv,data){
#    pred = sol.binomial.mixt(a, b, p.conv*data$p.factor[1], data$error[1], data)
    pred = sol.binomial.mixt(a, b, p.conv, data$error[1], data)
    res = data$akn
    sum(sum(data$akn) * (pred - res)^2)
}

cost.function.per.T.error <- function (params, data, p.error) 
{
    a = params[1]
    b = 1 - a
    p.conv = params[2]
    ll <- lapply(data,function(x){binomial.solution.wrapper(a,b,p.conv,x)})
    norm <- sum(unlist(lapply(data,function(x){x$sum.akn[1]})))
    sum(unlist(ll),na.rm=T)/norm
}



call.opt.per.T.error <- function (params, data, data.error) 
{
    data <- data %>% group_by(N) %>% mutate(sum.akn = sum(akn))
    data %>%
        group_by(N) %>%
        do(ll=as_data_frame(.[c("N","k","akn","error","sum.akn","p.factor")])) %$%
    setNames(ll,N) ->
    data

    oo <- optim(par = params, fn = cost.function.per.T.error, 
        lower = c(10^(-6), 0.001), upper = c(1, 0.5), method = "L-BFGS-B", 
        data = data)
    oo$par <- c(oo$par)
    names(oo$par) <- c("a", "p.conv")
    as_data_frame(t(as.matrix(c(oo$par, cost = oo$value))))
}


generate.random.per.T.error <- function (params, data, no.random, random.fact = 0.5, no.cores = 1) 
{
    rnd.lhs = random.fact + (1/random.fact - random.fact) * lhs::randomLHS(n = no.random, 
        k = length(params))
    df.0 <- t(apply(rnd.lhs, 1, function(x) {
        x * params
    }))
    opt.out <- mclapply(1:nrow(df.0), function(x) {
        call.opt.per.T.error(params = df.0[x, ], data = data)
    }, mc.cores = no.cores)
    bind_rows(opt.out, .id = "run") %>% top_n(1, -cost)
}


multistep.per.T.error <- function (data, data.error, threshold.excluded.events = 25,  no.random, random.fact) 
{
    data.error %>%
        filter(N>0) %>%
        group_by(N) %>%
        do(error=ml.estimate.pconv(data=.)) %>%
        unnest ->
        data.error
    p.error=median(data.error$error,na.rm=T)
    all.conversions <- compute.conversion.rate(data, p.error = p.error, 
        threshold.excluded.events = threshold.excluded.events)
    p.conv = ifelse(is.na(all.conversions$pr), 0.5, all.conversions$pr)
    mixture.fit.first <- call.opt.mixture(params = 0.5, data = data, 
                                          p.error = p.error, p.conv = p.conv)
    colnames(data)=c("N","k","akn")
    data.full=left_join(data,data.error,by="N")
    data.full %>%
        ungroup %>%
        mutate(p.error.min=quantile(error,0.1,na.rm=T)) %>%
        mutate(p.factor=error/p.error.min)->
        data.full
    mixture.fit <- generate.random.per.T.error(c(mixture.fit.first$a, 
        p.conv), data = data.full, no.random = no.random, 
        random.fact = random.fact)
}



call.opt.mixture <- function(params,p.conv,p.error,data){
#    data <- data[,c("X1","X2","X4")]
    colnames(data) <- c("N","k","akn")
    data %>%
        group_by(N) %>%
        mutate(sum.akn=sum(akn))->
        data
    oo <- optim(par=params,fn=cost.fn.binomial,lower=c(10^(-6)),upper=c(1),method="L-BFGS-B",data=data,p.conv=p.conv,p.error=p.error)
    oo$par <- c(oo$par,1-oo$par)
    names(oo$par) <- c("a","b")
    as_data_frame(t(as.matrix(c(oo$par,cost=oo$value))))
}


provide.total.number.of.events.per.N <- function(data,lower.N=20,upper.N=50){
    colnames(data) <- c("X1","X2","X3")
    data %>%
        filter(X1>lower.N & X1< upper.N) %>%
        group_by(X1) %>%
        summarise(no.events=sum(X3)) %>%
        ungroup %>%
        filter(no.events>0)
}
