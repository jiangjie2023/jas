

################################################################################
# DISTRIBUTION FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION (BASED ON library normalp).
pnormp=function (q, mu = 0, sigmap = 1, p = 2, lower.tail = TRUE, log.pr = FALSE) 
{
    if (!is.numeric(q) || !is.numeric(mu) || !is.numeric(sigmap) || 
        !is.numeric(p)) 
        stop(" Non-numeric argument to mathematical function")
    if (min(p) < 1) 
        stop("p must be at least equal to one")
    if (min(sigmap) <= 0) 
        stop("sigmap must be positive")
    z <- (q - mu)/sigmap

    zz <- abs(z)^p
    zp <- pgamma(zz, shape = 1/p, scale = p)
    lzp <- pgamma(zz, shape = 1/p, scale = p, log=TRUE)
    zp <- ifelse(z < 0, 0.5 - exp(lzp-log(2)), 0.5 + exp(lzp-log(2)))
    if (log.pr == TRUE) 
    zp<-ifelse(z < 0,log(0.5-exp(lzp-log(2))), log(0.5+exp(lzp-log(2))))
    zp
}

################################################################################
# DENSITY FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION (BASED ON library normalp).
dnormp=function(x, mu = 0, sigmap = 1, p = 2, log = FALSE) 
{
    if (!is.numeric(x) || !is.numeric(mu) || !is.numeric(sigmap) || 
        !is.numeric(p)) 
        stop(" Non-numeric argument to mathematical function")
    if (min(p) < 1) 
        stop("p must be at least equal to one")
    if (min(sigmap) <= 0) 
        stop("sigmap must be positive")
    cost <- 2 * p^(1/p) * gamma(1 + 1/p) * sigmap
    expon1 <- (abs(x - mu))^p
    expon2 <- p * sigmap^p
    dsty <- (1/cost) * exp(-expon1/expon2)
    if (log == TRUE) 
        dsty <- log(dsty)
    dsty
}

################################################################################
# LOG-LIKELIHOOD FUNCTION (REQUIRED FOR SEVERAL .LEP FUNCTIONS)
log.lik.LEP<-function(y,Cens,mu,sigma,alpha)
{
	n=length(y); aux<-rep(0,n) 
	SP=as.vector(sigma*(1/alpha)^(1/alpha))	
	 
aux=Cens*(dnormp(log(y),mu,sigmap=SP,p=alpha,log=TRUE)-log(y))+(1-Cens)*log(1-pnormp(log(y),mu,sigmap=SP,p=alpha))
	 
	return(sum(aux))
}




################################################################################
# DIC
DIC.LEP=function(y,Cens, chain)
{
	chain=as.matrix(chain); N=dim(chain)[1];n=length(y); LL<-rep(0,times=N)

	for(iter in 1:N)
	{
LL[iter]=log.lik.LEP(y,Cens, mu= chain[iter,1],
sigma=chain[iter,2],alpha=chain[iter,3])
	}

	aux=apply(chain,2,"mean")
	pd=-2*mean(LL)+2*log.lik.LEP(y,Cens, mu= chain[iter,1],
sigma=chain[iter,2],alpha=chain[iter,3])
 
	DIC=-2*mean(LL)+pd

	 
	return(DIC)  
}
#DIC.LEP(y,Cens,MCMCsamples)