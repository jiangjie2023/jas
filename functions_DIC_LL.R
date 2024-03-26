


library(VGAM)


################################################################################
# LOG-LIKELIHOOD FUNCTION (REQUIRED FOR SEVERAL .LLAP FUNCTIONS)
log.lik.LLAP<-function(y,Cens,mu,sigma)
{
  n=length(y); aux<-rep(0,n)
  
  aux=Cens*(dlaplace(log(y),mu,sigma,log=TRUE)-log(y))+
    (1-Cens)*log(1-plaplace(log(y),mu, sigma))
  
  return(sum(aux))
}
 

################################################################################
# DIC
DIC.LLAP=function(y,Cens, chain)
{
  chain=as.matrix(chain); N=dim(chain)[1];n=length(y); LL<-rep(0,times=N)
  
  for(iter in 1:N)
  {
    LL[iter]=log.lik.LLAP(y,Cens, mu= chain[iter,1],sigma=chain[iter,2])
  }
  
  aux=apply(chain,2,"mean")
  pd=-2*mean(LL)+2*log.lik.LLAP(y,Cens, mu= chain[iter,1],sigma=chain[iter,2])
  
  DIC=-2*mean(LL)+pd
  return(DIC)  
}
