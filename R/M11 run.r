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
           - theta)*(y<=200)
           +(normal_lpdf(y|theta,sigob))*(y>200);

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
real b8;
real b9;
real b10;
real<lower=0> sig2;
vector[N] e;
real sigob;
}


#### define some parameters which is the function of other parameters 
transformed parameters {
vector[N] mu;

### define the mean function


mu[1] <- y[1];  ## not tested

  for (n in 2:N) {
	
      mu[n] <- mu[n-1]*exp(u+b8*((1-exp(-n/b10))/(n/b10))+b9*((1-exp(-n/b10))/(n/b10)-exp(-n/b10))+e[n]); 
	    
	  
	  }

}



##### define piror and likelihood function
model {
lambda~uniform(-1,1);
u ~ normal(0.028,0.01);
b8 ~ normal(0.05,0.01);
b9 ~ normal(0.05,0.01);
b10 ~ normal(0.05,0.01);
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

##  y=read.csv("COVID-19_confirmed__Germany.csv")
##  y=read.csv("COVID-19_confirmed__Italy.csv")
## y=read.csv("COVID-19_confirmed__Japan.csv")
##  y=read.csv("COVID-19_confirmed__Spain.csv")
##  y=read.csv("COVID-19_confirmed__United Kingdom.csv")
  y=read.csv("COVID-19_confirmed__US.csv")
##  y=read.csv("COVID-19_confirmed_Australia.csv")

y=y[,2]




E=rep(1,length(y))

N=length(y)
dat <- list(N = N, y = y,E=E)


#### initial values
inits =list(list(lambda=0,u=0.028,b8=0.05,b9=0.05,b10=0.05,sig2=0.0005,e=rep(0.01,N),sigob=1000))  ### 134


fit <- stan(model_code = stanmodelcode, model_name = "example",init =inits,
data = dat, warmup=10000,iter = 100000, chains = 1, seed=999,thin=1)
print(fit,digits = 4)














####################  trace plot

library(methods)
library(rstan)


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")


#k=7

pdf(paste("M11_",abn[k],"_p1.pdf",sep=""), width = 16.5, height = 8.50) 

traceplot(fit, pars=c("u","b8","b9","b10"), inc_warmup = T,  nrow = 4, ncol = 1, window = NULL, include = TRUE)

dev.off()

pdf(paste("M11_",abn[k],"_p2.pdf",sep=""), width = 16.5, height = 8.50) 

traceplot(fit, pars=c("sigob","sig2","lambda"), inc_warmup = T,  nrow = 3, ncol = 1, window = NULL, include = TRUE)

dev.off()



##################   plot 30/45/60/   length(y)


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , mean )

Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , quantile, probs=0.05 )

Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , quantile, probs=0.95 )






pdf(paste("M11_",abn[k],"_30.pdf",sep=""), width = 16.5, height = 8.50) 

m=30

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M11_",abn[k],"_45.pdf",sep=""), width = 16.5, height = 8.50) 

m=45

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M11_",abn[k],"_60.pdf",sep=""), width = 16.5, height = 8.50) 

m=60

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M11_",abn[k],"_all.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),lwd = 2)

dev.off()






########    log plot          ####################################################

k=1


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")



Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , mean )

Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , quantile, probs=0.05 )

Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , quantile, probs=0.95 )




y=log(y)
Epars_m=log(Epars_m)
Epars_9=log(Epars_9)
Epars_5=log(Epars_5)


setwd("./covid19modelrisk/output/plot")

pdf(paste("M11_",abn[k],"_log.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M11 model for ",name[k],sep=""),lwd = 2)

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

################################################ as.data.frame(fit@sim[[1]][[1]])[1,]

P=20
Mu=matrix(,nrow=90000,ncol=P)
X=matrix(,nrow=90000,ncol=P)


lambda_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,1]
u_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,2]
b8_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,3]
b9_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,4]
b10_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,5]
sig2_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,6]


set.seed(999)

for (r in 1:90000){



Mu[r,1]=y[length(y)]
for (s in 1:(P-1) ){

e=rnorm(1,0,sig2_p[r])

Mu[r,s+1]=Mu[r,s]*exp(u_p[r]+b8_p[r]*((1-exp(-Mu[r,s]/b10_p[r]))/(Mu[r,s]/b10_p[r]))+b9_p[r]*((1-exp(-Mu[r,s]/b10_p[r]))/(Mu[r,s]/b10_p[r])-exp(-Mu[r,s]/b10_p[r]))+e)



}}




for (r in 1:90000){
for (s in 1:(P-1) ){

X[r,s]=GP.r(Mu[r,s]*(1-lambda_p[r]), lambda_p[r])


}}


#X=Mu


Xm=apply(X,2,mean)

X5 <- apply( X , 2 , quantile, probs=0.05 )

X9 <- apply( X , 2 , quantile, probs=0.95 )



Zm=c(y,Xm)
Z5=c(y,X5)
Z9=c(y,X9)


k=1

name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")


m=length(Zm)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(Zm[1:m],ylim=c(min(Zm[1:m]),max(Zm[1:m],Z9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Z9[1:m]), Z5[1:m]), col = "grey80", border = NA)

lines(newx,Z9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Z5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(Zm[1:m],ylim=c(min(Zm[1:m]),max(Zm[1:m],Z9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19)














########################################################################################################################

n=length(y)



u=0.028
b8=0.05
b9=0.05
b10=0.05


x=c()
x[1]=100
for (i in 1:(n-1) ){

x[i+1]=x[i]*exp(u+b8*((1-exp(-x[i]/b10))/(x[i]/b10))+b9*((1-exp(-x[i]/b10))/(x[i]/b10)-exp(-x[i]/b10)))
}
x=round(x)
x



plot(y,ylim=c(min(y,x),max(y,x)))
par(new=TRUE)
plot(x,col="red",ylim=c(min(y,x),max(y,x)))









