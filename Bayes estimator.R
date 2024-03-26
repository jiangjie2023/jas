########改变a,A,B,A0,B0的值对应不同的先验分布

library(truncdist)
library(truncnorm)
library(DIRECT)
library(digest)
library(miscTools)
library(zoo)

mu=2
sigma=0.5
alpha=3
  
censor=10
a=1
##MH, Using Gaussian random walk proposals
numb=100
N=c(100,150,200,300)

se=sd_0=be_0=mse_be=bias_be=array(0,c(length(N),3))
BE=array(0,c(3,numb,length(N)))
cp_mu=cp_sigma=cp_alpha=c(0)
for(ii in 1:length(N))
{

n=N[ii]


cp1=cp2=cp3=0

 be=array(0,c(3,numb)) 

for(k in 1:numb)
{

logy=c(0)
for(i in 1:n)
{
  u=rgamma(1,1+1/alpha,1)
  logy[i]=runif(1,mu-sigma*u^(1/alpha), mu+sigma*u^(1/alpha))
}
y=exp(logy)
cen_ind=sample(1:n,(censor/100*n))
y_data=y[cen_ind]
C_data=c(0)
for(i in 1:length(cen_ind))
{C_data[i]=runif(1,0,y_data[i])
}
C=y
C[cen_ind]=C_data#????????min(y,C)
Cens=rep(1,length(y))
Cens[cen_ind]=0

mu0=2
sigma0=0.5
alpha0=3
stepsize=c(1,1,0.015)
tuning=50
burn_in=10000
burn=matrix(0,3,burn_in)
for(i in 1:burn_in)
{ 
u_curr=c(0)
for(j1 in 1:n){
u_curr[j1]=rtrunc(1,spec="exp",a=(abs(log(y[j1])-mu0)/sigma0)^alpha0,b=Inf,rate=1)
}
u0=u_curr

#logy_n0=c(0)
#for(i1 in 1:(censor/100*n))
#{ 
#b1=mu0+sigma0*u0[cen_ind[i1]]^(1/alpha0)
#if(log(C_data[i1])>b1){logy_n0[i1]=runif(1,mu0-sigma0*u0[cen_ind[i1]]^(1/alpha0), mu0+sigma0*u0[cen_ind[i1]]^(1/alpha0))}else
 # {
#logy_n0[i1]=rtrunc(1,spec="unif",a=log(C_data[i1]),b=Inf,mu0-sigma0*u0[cen_ind[i1]]^(1/alpha0), mu0+sigma0*u0[cen_ind[i1]]^(1/alpha0))  
#}
#}
#y[cen_ind]<-exp(logy_n0)

a_1=apply(cbind(mu0-sigma0*u0^(1/alpha0),log(C)),1,max)
b_1=as.vector(mu0+sigma0*u0^(1/alpha0))
logy=Cens*log(C)+
(1-Cens)*(as.numeric(I(a_1<b_1))*runif(n,min=apply(cbind(a_1,b_1),1,min),max=apply(cbind(a_1,b_1),1,max))+(1-as.numeric(I(a_1<b_1)))*log(C))	
y<-exp(logy)

  propmu=rnorm(1,mu0,stepsize[1])
  a0=exp(-1/sigma0^alpha0*sum(abs(log(y)-propmu)^alpha0))/exp(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0)
  ) 
  if(runif(1)<a0){mu0=propmu}else{mu0=mu0}
 burn[1,i]=mu0
  props=rtruncnorm(1,0.001,Inf,sigma0,stepsize[2])
  #a1=props^(-n-a)*exp(-1/props^alpha0*sum(abs(log(y)-mu0)^alpha0))/(sigma0^(-n-a)*exp(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0)))
 a1= (-n-a)*log(props)-1/props^alpha0*sum(abs(log(y)-mu0)^alpha0)-
((-n-a)*log(sigma0)-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0))
if(a1< -100){a1=-100}
if(runif(1)<exp(a1)){sigma0=props}else{sigma0=sigma0}
 burn[2,i]=sigma0
  A=1+digamma(1+1/alpha0)
  B=(1+1/alpha0)*trigamma(1+1/alpha0)-1
  propalpha=rtruncnorm(1,1,Inf,alpha0,stepsize[3])
  A0=1+digamma(1+1/propalpha)
  B0=(1+1/propalpha)*trigamma(1+1/propalpha)-1
  a2= 1/gamma(1+1/propalpha)^n*exp(-1/sigma0^propalpha*sum(abs(log(y)-mu0)^propalpha)
  )*((A0^2+B0)*propalpha^(-3))^(1/2)/(1/gamma(1+1/alpha0)^n*exp(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0)
  )*((A^2+B)*alpha0^(-3))^(1/2))
  if(runif(1)<a2)
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

num=10000
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
  a0=exp(-1/sigma0^alpha0*sum(abs(log(y)-propmu)^alpha0))/exp(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0)
  ) 
  if(runif(1)<a0){mu0=propmu;ac0=ac0+1}else{mu0=mu0}
  
  props=rtruncnorm(1,0,Inf,sigma0,stepsize[2])
 a1= (-n-a)*log(props)-1/props^alpha0*sum(abs(log(y)-mu0)^alpha0)-
((-n-a)*log(sigma0)-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0))
if(a1< -100){a1=-100}
if(runif(1)<exp(a1)){sigma0=props;ac1=ac1+1}else{sigma0=sigma0}
 A=1+digamma(1+1/alpha0)
  B=(1+1/alpha0)*trigamma(1+1/alpha0)-1
  propalpha=rtruncnorm(1,1,Inf,alpha0,stepsize[3])
  A0=1+digamma(1+1/propalpha)
  B0=(1+1/propalpha)*trigamma(1+1/propalpha)-1
  a2= 1/(gamma(1+1/propalpha))^n*exp(-1/sigma0^propalpha*sum(abs(log(y)-mu0)^propalpha)
  )*((A0^2+B0)*propalpha^(-3))^(1/2)/(1/(gamma(1+1/alpha0))^n*exp(-1/sigma0^alpha0*sum(abs(log(y)-mu0)^alpha0)
  )*((A^2+B)*alpha0^(-3))^(1/2))
  if(runif(1)<a2)
{alpha0=propalpha;ac2=ac2+1}else{alpha0=alpha0}
  Mchian[j,]=c(mu0,sigma0,alpha0)
#if(j%%1000==0){print(j)}
}
MCMCsamples=Mchian[seq(1,num,by=50),]
be[,k]=apply(MCMCsamples, 2,mean)#######################?玫?be

palpha_mu=quantile(MCMCsamples[,1],probs=0.95)
if(mu<palpha_mu){cp1=cp1+1}
palpha_sigma=quantile(MCMCsamples[,2],probs=0.95)
if(sigma<palpha_sigma){cp2=cp2+1}

palpha_alpha=quantile(MCMCsamples[,3],probs=0.95)
if(alpha<palpha_alpha){cp3=cp3+1}

}
cp_mu[ii]=cp1
cp_sigma[ii]=cp2
cp_alpha[ii]=cp3
BE[,,ii]=be
be_0[ii,]=apply(be,1,mean)
sd_0[ii,]=apply(be,1,sd)
se[ii,]=sd_0[ii,]/sqrt(numb)

real=c(mu,sigma,alpha)
 MSE_be=c(0)
for(i in 1:3)
{
MSE_be[i]=mean((be[i,]-real[i])^2)
}
mse_be[ii,]=MSE_be
 
Bias_be=c(0)
for(i in 1:3)
{
Bias_be[i]=mean(be[i,]-real[i])
}
bias_be[ii,]=Bias_be

}


be_0
se
mse_be
bias_be




##########################################
########################################
###################################mse plot####################
 
library(tidyverse)
library(ggplot2)
library(ggpubr)
N=c(100,150,200,300)

r1=mse1[,2];r2=mse2[,2];r3=mse3[,2];r4=mse4[,2];r5=mse5[,2]
Prior=rep(c("R1","R2","R3","R4","RT"),each=4)
df<-data.frame(MSE=c(r1,r2,r3,r4,r5),Prior,n=N)
ggplot(data=df,mapping = aes(x =n, y = MSE,linetype=Prior,
 colour=Prior,shape=Prior)) + geom_line(size=0.8)+ geom_point(size=2)+
theme_bw() +#去掉背景灰色
theme(plot.title=element_text(hjust=0.5))+#让标题题目居中
labs(title=expression(sigma))+
 theme(
        legend.position = c(0.85,0.70))+#更改图例的位置，放至图内部的右上角
scale_color_manual(values= c("#007896FF","#619CFF","#3AC96DFF","#ef1828","#00B0F6"))+ #自定义点颜色
  scale_shape_manual(values= c(15,17,19,21,23)) #自定义形状

 
 ##########################################
########################################
###################################bias plot###################

r1=bias1[,2];r2=bias2[,2];r3=bias3[,2];r4=bias4[,2];r5=bias5[,2]
Prior=rep(c("R1","R2","R3","R4","RT"),each=4)
df<-data.frame(Bias=c(r1,r2,r3,r4,r5),Prior,n=N)
ggplot(data=df,mapping = aes(x =n, y = Bias,linetype=Prior,
 colour=Prior,shape=Prior)) + geom_line(size=0.8)+ geom_point(size=2)+ylim(-0.2,0.2)+
geom_hline(yintercept=0.0,color="purple")+
theme_bw() +#去掉背景灰色
theme(plot.title=element_text(hjust=0.5))+#让标题题目居中
labs(title=expression(sigma))+
 theme(
        legend.position = c(0.85,0.75))+#更改图例的位置，放至图内部的右上角
scale_color_manual(values= c("#007896FF","#619CFF","#3AC96DFF","#ef1828","#00B0F6"))+ #自定义点颜色
  scale_shape_manual(values= c(15,17,19,21,23)) #自定义点颜色

 