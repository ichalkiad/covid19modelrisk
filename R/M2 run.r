###############################################################################
# "Infection Rate Models for COVID-19: 
#  Model Risk and Public Health News Sentiment Exposure Adjustments"
#   
#  Ioannis Chalkiadakis, Kevin Hongxuan Yan, Gareth W. Peters, Pavel V. Shevchenko
#
#  Kevin Hongxuan Yan
#  March 2021
###############################################################################


rm(list=ls(all=TRUE))



library(methods)
library(rstan)


cat("Stan version:", stan_version(), "\n")
stanmodelcode <- "


#### define GP function
functions {

real normal_generalized_poisson_log(int y, real theta, real lambda, real sigob) { 
    return (log(theta) 
           + (y - 1) * log(theta + y * lambda) 
           - lgamma(y + 1) 
           - y * lambda 
           - theta)*(y<=500)
           +(normal_lpdf(y|theta,sigob))*(y>500);
		  ### +log(1-normal_cdf(0.08*u,u,0.01))*(b0>0.08*u);   #### +log(b0>0.08*u); (b0>0.08*u)
  } 

}


#### define our data, integer value
data {
int<lower=0> N;
int y[N];
real E[N];
}


#### define parameters  
parameters {

real lambda;
real u;
real b0;
real<lower=0> sig2;
vector[N] e;
real sigob;
}


#### define some parameters which is the function of other parameters 
transformed parameters {
vector[N] mu;

### define the mean function
#mu[1] <- exp(u*exp(-b0)+e[1]);

mu[1] <- y[1];  ## not tested

  for (n in 2:N) {
	
      mu[n] <- mu[n-1]*exp(u*exp(-b0*n)+e[n]); 
	  
	  }

}



##### define piror and likelihood function
model {
lambda~uniform(-1,1);
u ~ normal(0.24,0.01);
b0 ~ normal(0.031,0.001);
sig2 ~  gamma(1,1); 
e~normal(0,sig2);
sigob ~ uniform(100,1100);

for (j in 1:N) {

y[j] ~ normal_generalized_poisson( mu[j]*E[j], lambda, sigob);
}


}



#### calculate DIC
generated quantities {
vector[N] log_lik;
real dev; // deviance


  for (n in 1:N){
    log_lik[n] <- normal_generalized_poisson_log (y[n],mu[n]*E[n], lambda, sigob);
  }
 
  dev <- 0;
        for ( n in 1:N ) {
            dev <- dev + (-2) * normal_generalized_poisson_log (y[n],mu[n]*E[n], lambda, sigob);
        }
  
}





"


###### 





####### read new data Y
setwd("./covid19modelrisk/data/covid19_infected/")

  y=read.csv("COVID-19_confirmed__Germany.csv")
##  y=read.csv("COVID-19_confirmed__Italy.csv")
##  y=read.csv("COVID-19_confirmed__Japan.csv")
##  y=read.csv("COVID-19_confirmed__Spain.csv")
##  y=read.csv("COVID-19_confirmed__United Kingdom.csv")
##  y=read.csv("COVID-19_confirmed__US.csv")
##  y=read.csv("COVID-19_confirmed_Australia.csv")

y=y[,2]




E=rep(1,length(y))

N=length(y)
dat <- list(N = N, y = y,E=E)


#### initial values
inits =list(list(lambda=0,u=0.24,b0=0.031,sig2=0.0005,e=rep(0.01,N),sigob=1000))  ### 134


fit <- stan(model_code = stanmodelcode, model_name = "example",init =inits,
data = dat, warmup=10000,iter = 100000, chains = 1, seed=999,thin=1)
print(fit,digits = 4)


####  threshold =500









###############################################################################







####################  trace plot

library(methods)
library(rstan)


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")




k=7


pdf(paste("M2_",abn[k],"_p1.pdf",sep=""), width = 16.5, height = 8.50) 

traceplot(fit, pars=c("u","b0","lambda"), inc_warmup = T,  nrow = 3, ncol = 1, window = NULL, include = TRUE)

dev.off()

pdf(paste("M2_",abn[k],"_p2.pdf",sep=""), width = 16.5, height = 8.50) 

traceplot(fit, pars=c("sigob","sig2"), inc_warmup = T,  nrow = 2, ncol = 1, window = NULL, include = TRUE)
dev.off()



##################   plot 30/45/60/   length(y)


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , mean )

Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , quantile, probs=0.05 )

Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , quantile, probs=0.95 )





pdf(paste("M2_",abn[k],"_30.pdf",sep=""), width = 16.5, height = 8.50) 

m=30

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M2_",abn[k],"_45.pdf",sep=""), width = 16.5, height = 8.50) 

m=45

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M2_",abn[k],"_60.pdf",sep=""), width = 16.5, height = 8.50) 

m=60

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M2_",abn[k],"_all.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),lwd = 2)

dev.off()



########    log plot          ####################################################

k=2


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")




Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , mean )

Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , quantile, probs=0.05 )

Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , quantile, probs=0.95 )


y=log(y)
Epars_m=log(Epars_m)
Epars_9=log(Epars_9)
Epars_5=log(Epars_5)

setwd("./covid19modelrisk/output/plot/")

pdf(paste("M2_",abn[k],"_log.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),lwd = 2)

dev.off()





#####################################################################################################################################################
######################################      forecast          ###################################

GP.r <- function(L, s)
{

v=0
Lp=exp(log(L)+ (v - 1) * log(L + v * s) - v * s - L - lgamma(v + 1))
sLp=Lp

# calculate the vector of cum pmf
v=1
while((sLp<0.9999999) & ((L + v * s )> 0))   
{Lp1=exp(log(L)+ (v - 1) * log(L + v * s) - v * s - L - lgamma(v + 1))
Lp=c(Lp,Lp1)
sLp=sLp+Lp1
v=v+1
}

csLp=c(cumsum(Lp),1) ##### add a max for it 
x=length(Lp)

# first check if you are in the first interval

U = runif(1)
B = FALSE
i = 1
      
while(B == FALSE) {    
    
  if(U > csLp[i])  ### U may greater than max of csLp if do not set 1
    { i=i+1}
  else
    {X=i-1
     B=TRUE}     

}

return(X)

}

################################################

P=20
Mu=matrix(,nrow=90000,ncol=P)
X=matrix(,nrow=90000,ncol=P)


lambda_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,1]
u_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,2]
b0_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,3]
sig2_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,4]


set.seed(999)

for (r in 1:90000){



Mu[r,1]=y[length(y)]
for (s in 1:(P-1) ){

e=rnorm(1,0,sig2_p[r])

Mu[r,s+1]=Mu[r,s]*exp(u_p[r]*exp(-b0_p[r]*(s+length(y)))+e)


}}



sigob=1000
for (r in 1:90000){
for (s in 1:P ){

X[r,s]=rnorm(1,Mu[r,s], sigob)


}}




Xm=apply(X,2,mean)

X5 <- apply( X , 2 , quantile, probs=0.05 )

X9 <- apply( X , 2 , quantile, probs=0.95 )



Zm=c(y,Xm)
Z5=c(y,X5)
Z9=c(y,X9)


##### read forecast
setwd("./covid19modelrisk/data/covid19_infected/")

##  FT=read.csv("COVID-19_confirmed__Germany_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed__Italy_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed__Japan_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed__Spain_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed__United Kingdom_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed__US_FT.csv")[,2]
FT=read.csv("COVID-19_confirmed_Australia_FT.csv")[,2]





M=length(FT)

e=FT-Xm
mae=mean(abs(e))
rmse=sqrt(mean(e^2))
############ PR
 p=e/FT
# p1=e1/mean(mo1)
mape=mean(abs(p))
########## QR
md=sum(abs(FT[2:M]-FT[1:(M-1)]))/(M-1)
qj=e/md
mase=mean(abs(qj))


list(mae,rmse,mape,mase)




Z_FT=c(y,FT)


k=7

name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")


m=length(Zm)


setwd("./covid19modelrisk/output/forecast/")

pdf(paste("M2_",abn[k],"_F.pdf",sep=""), width = 16.5, height = 8.50) 

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(Zm[1:m],ylim=c(min(Zm[1:m]),max(Zm[1:m],Z9[1:m])),ylab="Count",xlab="Time",main=paste("Out-sample forecast results of M2 model for ",name[k],sep=""),pch=19,col="red",xaxt = 'n')
polygon(c(rev(newx), newx), c(rev(Z9[1:m]), Z5[1:m]), col = "grey80", border = NA)

par(new=TRUE)
plot(Zm[1:m],ylim=c(min(Zm[1:m]),max(Zm[1:m],Z9[1:m])),ylab="Count",xlab="Time",main=paste("Out-sample forecast results of M2 model for ",name[k],sep=""),pch=19,col="red",xaxt = 'n')
lines(newx,Z9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Z5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(Z_FT[1:m],ylim=c(min(Zm[1:m]),max(Zm[1:m],Z9[1:m])),ylab="Count",xlab="Time",main=paste("Out-sample forecast results of M2 model for ",name[k],sep=""),pch=19,xaxt = 'n')
legend("topleft", legend=c("Observed data", "Forecasted data"), col=c("black", "red"), lty=c(1,1), lwd = 3,cex=2)

# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/27","Mar/17","May/05","Jun/25","Aug/06"))   ### GM
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/31","Mar/13","May/09","Jun/29","Aug/06"))   ### IT
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/22","Mar/12","Apr/30","Jun/19","Aug/06"))   ### JP
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Feb/01","Mar/20","May/10","Jun/30","Aug/06"))   ### SP
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/31","Mar/21","May/09","Jun/29","Aug/06"))   ### UK
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/22","Mar/12","Apr/30","Jun/19","Aug/06"))   ### US
axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/26","Mar/16","May/04","Jun/24","Aug/06"))   ### AU

dev.off()









######   RMSE results        ###############################################################

### Germany
> list(mae,rmse,mape,mase)
[[1]]
[1] 12517.73

[[2]]
[1] 14479.24

[[3]]
[1] 0.05425658

[[4]]
[1] 10.42504


## Italy
> list(mae,rmse,mape,mase)
[[1]]
[1] 28633.06

[[2]]
[1] 33356.15

[[3]]
[1] 0.1114455

[[4]]
[1] 42.5554


## Japan
> list(mae,rmse,mape,mase)
[[1]]
[1] 3407.51

[[2]]
[1] 3475.616

[[3]]
[1] 0.06044399

[[4]]
[1] 3.404822



## Spain
> list(mae,rmse,mape,mase)
[[1]]
[1] 6404.752

[[2]]
[1] 7764.212

[[3]]
[1] 0.01816862

[[4]]
[1] 1.153605



## UK
> list(mae,rmse,mape,mase)
[[1]]
[1] 21333.73

[[2]]
[1] 24918.48

[[3]]
[1] 0.06618652

[[4]]
[1] 20.42946



## US
> list(mae,rmse,mape,mase)
[[1]]
[1] 378799.5

[[2]]
[1] 464800.3

[[3]]
[1] 0.06784786

[[4]]
[1] 8.161634


## AU
> list(mae,rmse,mape,mase)
[[1]]
[1] 3069.551

[[2]]
[1] 3986.055

[[3]]
[1] 0.1254013

[[4]]
[1] 12.61277






##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
n=length(y)

u=0.04
b0=0.004



n=length(y)

x=c()
x[1]=100
for (i in 1:(n-1) ){

x[i+1]=x[i]*exp(u*exp(-b0*i))
}
x=round(x)


plot(y,ylim=c(min(y),max(y)))
par(new=TRUE)
plot(x,col="red",ylim=c(min(y),max(y)))




#####   initial values         ##########################################################################

### Germany
u=0.24
b0=0.031#u*0.08 ### 


## Italy
u=0.24
b0=0.03#u*0.08 ### 



## Japan
u=0.05
b0=0.0055



## Spain
u=0.3
b0=0.037



## UK
u=0.3
b0=0.037



## US
u=0.1
b0=0.0067


## AU
u=0.04
b0=0.004



















