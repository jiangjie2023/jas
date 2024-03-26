
library(truncdist)
library(truncnorm)
library(DIRECT)
library(maxLik) 
library(rethinking)
y=scan()
4.50 19.13 14.24  7.87  5.49  2.02  9.22  3.82 26.31  4.65  2.62  0.90
21.73  0.87  0.51  3.36 43.01  0.81  3.36  1.46 24.80 10.86 17.14 15.96
7.28  4.33 22.69  2.46  3.48  4.23  6.54  8.65   5.41 2.23 4.34 32.15  4.87  5.71  7.59
3.02  4.51  1.05  9.47 79.05  2.02  4.26 11.25 10.34 10.66 12.03  2.64
14.76  1.19  8.66 14.83  5.62 18.10 25.74 17.36  1.35  9.02  6.94  7.26
4.70  3.70  3.64  3.57 11.64 6.25 25.82  3.88  3.02 19.36 20.28 46.12  5.17  0.20 36.66
10.06  4.98  5.06 16.62 12.07  6.97  0.08  1.40  2.75  7.32  1.26  6.76
8.60  7.62  3.52  9.74  0.40  5.41  2.54  2.69  8.26  0.50  5.32  5.09 2.09 7.93 12.02
13.80  5.85  7.09  5.32  4.33  2.83  8.37 14.77  8.53 11.98  1.76  4.40
34.26  2.07 17.12 12.63  7.66  4.18 13.29 23.63  3.25  7.63  2.87  3.31
2.26  2.69 11.79  5.34  6.93 10.75 13.11  7.39

cen_ind=c(10,14,21,22,64,72,73,91,110)
C=y  
Cens=rep(1,length(y))
Cens[cen_ind]=0

a=1#与先验有关
n=137
mu0=4
sigma0=1.5
alpha0=2
stepsize=rep(1,3)
tuning=1000
burn_in=100000
mu1=c(0)
burn=matrix(0,3,burn_in)
for(i in 1:burn_in)
{ 
  u_curr=c(0)
  for(j1 in 1:n){
    u_curr[j1]=rtrunc(1,spec="exp",a=(abs(log(y[j1])-mu0)/sigma0)^alpha0,b=Inf,rate=1)
  }
  u0=u_curr
  
  a_1=apply(cbind(mu0-sigma0*u0^(1/alpha0),log(C)),1,max)
  b_1=as.vector(mu0+sigma0*u0^(1/alpha0))
  logy=Cens*log(C)+
    (1-Cens)*(as.numeric(I(a_1<b_1))*runif(n,min=apply(cbind(a_1,b_1),1,min),max=apply(cbind(a_1,b_1),1,max))+(1-as.numeric(I(a_1<b_1)))*log(C))	
  y<-exp(logy)
  
  
  propmu=rnorm(1,mu0,stepsize[1])
  mu1[i]=propmu
  #a0=exp(-1/sigma0^alpha0*sum(abs(log(y)-propmu)^alpha0))/exp(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0))
  a0= -1/sigma0^alpha0*sum(abs(log(y)-propmu)^alpha0)-(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0))
  
  if(a0< -100){a0=-100}
  if(runif(1)<exp(a0)){mu0=propmu}else{mu0=mu0}
  burn[1,i]=mu0
  props=rtruncnorm(1,0,Inf,sigma0,stepsize[2])
  #a1=props^(-n-a)*exp(-1/props^alpha0*sum(abs(log(y)-mu0)^alpha0)
  #)/(sigma0^(-n-a)*exp(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0) ))
  
  a1=(-n-a)*props-1/props^alpha0*sum(abs(log(y)-mu0)^alpha0)-
    ((-n-a)*sigma0-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0))
  if(a1< -100){a1=-100}
  if(runif(1)<exp(a1)){sigma0=props}else{sigma0=sigma0}
  burn[2,i]=sigma0
 A=1+digamma(1+1/alpha0)
  # A=0
  B=(1+1/alpha0)*trigamma(1+1/alpha0)-1
  propalpha=rtruncnorm(1,1,Inf,alpha0,stepsize[3])
  A0=1+digamma(1+1/propalpha)
  B0=(1+1/propalpha)*trigamma(1+1/propalpha)-1
  # a2= 1/gamma(1+1/propalpha)^n*exp(-1/sigma0^propalpha*sum(abs(log(y)-mu0)^propalpha)
  #)*((A0^2+B0)*propalpha^(-3))^(1/2)/(1/gamma(1+1/alpha0)^n*exp(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0)
  #)*((A^2+B)*alpha0^(-3))^(1/2))
  
  a2=-n*log(gamma(1+1/propalpha))-1/sigma0^propalpha*sum(abs(log(y)-mu0)^propalpha)+1/2*log((A0^2+B0)*propalpha^(-3))-
    (-n*log(gamma(1+1/alpha0))-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0)+1/2*log((A^2+B)*alpha0^(-3)))
  if(a2< -100){a2=-100}
  if(runif(1)<exp(a2))
  {alpha0=propalpha }else{alpha0=alpha0}
  burn[3,i]=alpha0
  if(i%%tuning==0){
    for(j in 1:3)
    {
      I0=i/tuning
      Index.range0=(tuning*(I0-1)+1):(tuning*I0-1)
      Index.range1=(tuning*(I0-1)+2):(tuning*I0)
      A.rate=1-mean(burn[j,Index.range0]==burn[j,Index.range1])
      if(A.rate>0.5){stepsize[j]=exp(log(stepsize[j])+0.01)}
      if(A.rate<0.3){stepsize[j]=exp(log(stepsize[j])-0.01)}
    }  
  }
}

num=100000
ac0=ac1=ac2=0
Mchian=array(0,c(num,3))
for(j in 1:num)
{
  u_curr=c(0)
  for(j1 in 1:n){
    u_curr[j1]=rtrunc(1,spec="exp",a=(abs(log(y[j1])-mu0)/sigma0)^alpha0,b=Inf,rate=1)
  }
  u0=u_curr
  
  a_1=apply(cbind(mu0-sigma0*u0^(1/alpha0),log(C)),1,max)
  b_1=as.vector(mu0+sigma0*u0^(1/alpha0))
  logy=Cens*log(C)+
    (1-Cens)*(as.numeric(I(a_1<b_1))*runif(n,min=apply(cbind(a_1,b_1),1,min),max=apply(cbind(a_1,b_1),1,max))+(1-as.numeric(I(a_1<b_1)))*log(C))	
  y<-exp(logy)
  
  
  propmu=rnorm(1,mu0,stepsize[1])
  #a0=exp(-1/sigma0^alpha0*sum(abs(log(y)-propmu)^alpha0))/exp(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0)) 
  a0= -1/sigma0^alpha0*sum(abs(log(y)-propmu)^alpha0)-(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0))
  if(a0< -100){a0=-100}
  
  if(runif(1)<exp(a0)){mu0=propmu;ac0=ac0+1}else{mu0=mu0}
  
  props=rtruncnorm(1,0,Inf,sigma0,stepsize[2])
  
  a1=(-n-a)*props-1/props^alpha0*sum(abs(log(y)-mu0)^alpha0)-
    ((-n-a)*sigma0-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0))
  if(a1< -100){a1=-100}
  
  if(runif(1)<exp(a1)){sigma0=props;ac1=ac1+1}else{sigma0=sigma0}
  A=1+digamma(1+1/alpha0)
  #A=0
  B=(1+1/alpha0)*trigamma(1+1/alpha0)-1
  propalpha=rtruncnorm(1,1,Inf,alpha0,stepsize[3])
  A0=1+digamma(1+1/propalpha)
  B0=(1+1/propalpha)*trigamma(1+1/propalpha)-1
  a2=-n*log(gamma(1+1/propalpha))-1/sigma0^propalpha*sum(abs(log(y)-mu0)^propalpha)+1/2*log((A0^2+B0)*propalpha^(-3))-
    (-n*log(gamma(1+1/alpha0))-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0)+1/2*log((A^2+B)*alpha0^(-3)))
  if(a2< -100){a2=-100}
  if(runif(1)<exp(a2))
  {alpha0=propalpha;ac2=ac2+1}else{alpha0=alpha0}
  Mchian[j,]=c(mu0,sigma0,alpha0)
  
  if(j%%1000==0){print(j)}
}



par(mfrow=c(3,1))
plot(Mchian[ ,1],ylab=expression(mu),type ="l",main="")
plot(Mchian[ ,2],ylab=expression(sigma),type ="l",main="")
plot(Mchian[,3],ylab=expression(alpha),type ="l",main="")


MCMCsamples=Mchian[seq(1,num,by=20),]
par(mfrow=c(3,1))
acf(MCMCsamples[,1],ylab=expression(mu))
acf(MCMCsamples[,2],ylab=expression(sigma))
acf(MCMCsamples[,3],ylab=expression(alpha))

be=apply(MCMCsamples, 2,mean)#######################得到be
sd=apply(MCMCsamples,2, sd)

y=scan()
4.50 19.13 14.24  7.87  5.49  2.02  9.22  3.82 26.31  4.65  2.62  0.90
21.73  0.87  0.51  3.36 43.01  0.81  3.36  1.46 24.80 10.86 17.14 15.96
7.28  4.33 22.69  2.46  3.48  4.23  6.54  8.65   5.41 2.23 4.34 32.15  4.87  5.71  7.59
3.02  4.51  1.05  9.47 79.05  2.02  4.26 11.25 10.34 10.66 12.03  2.64
14.76  1.19  8.66 14.83  5.62 18.10 25.74 17.36  1.35  9.02  6.94  7.26
4.70  3.70  3.64  3.57 11.64 6.25 25.82  3.88  3.02 19.36 20.28 46.12  5.17  0.20 36.66
10.06  4.98  5.06 16.62 12.07  6.97  0.08  1.40  2.75  7.32  1.26  6.76
8.60  7.62  3.52  9.74  0.40  5.41  2.54  2.69  8.26  0.50  5.32  5.09 2.09 7.93 12.02
13.80  5.85  7.09  5.32  4.33  2.83  8.37 14.77  8.53 11.98  1.76  4.40
34.26  2.07 17.12 12.63  7.66  4.18 13.29 23.63  3.25  7.63  2.87  3.31
2.26  2.69 11.79  5.34  6.93 10.75 13.11  7.39

DIC=DIC.LEP(y,Cens,MCMCsamples)

be
sd/sqrt(length(MCMCsamples))
DIC
HPDI(Mchian[,1], prob=0.95)
HPDI(Mchian[,2], prob=0.95)
HPDI(Mchian[,3], prob=0.95)
