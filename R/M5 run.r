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

real normal_generalized_poisson_log(int y, real theta, real lambda, real sigob, real u, real b1) { 
    return (log(theta) 
           + (y - 1) * log(theta + y * lambda) 
           - lgamma(y + 1) 
           - y * lambda 
           - theta)*(y<=200)
           +(normal_lpdf(y|theta,sigob))*(y>200)
		   +log(1-normal_cdf(u*-0.00000005,u,0.01))*(b1<u*-0.00000005);   #### +log(b0>0.08*u); (b0>0.08*u)
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
real b1;
real b3;
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
	
      mu[n] <- (mu[n-1]^2)*exp(u+b1*mu[n-1]+e[n])/(b3+mu[n-1]); 
 
 
 
	  }

}



##### define piror and likelihood function
model {
lambda~uniform(-1,1);
u ~ normal(0.1,0.01);
b1 ~ normal(-0.0000005,0.000000001);
b3 ~ normal(0,0.000001);
sig2 ~  gamma(1,1); 
e~normal(0,sig2);
sigob ~ uniform(100,1100);

for (j in 1:N) {

y[j] ~ normal_generalized_poisson( mu[j]*E[j], lambda, sigob,u,b1);
}


}


"


###### 





####### read new data Y
setwd("./covid19modelrisk/output/forecast")

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
inits =list(list(lambda=0,u=0.1,b1=-0.0000005,b3=0,sig2=0.01,e=rep(0.01,N),sigob=1000))  ### 134


fit <- stan(model_code = stanmodelcode, model_name = "example",init =inits,
data = dat, warmup=4000,iter = 40000, chains = 1, seed=999,thin=1)
print(fit,digits = 4)






Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[4001:40000,(6+N):(5+2*N)] , 2 , mean )

Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[4001:40000,(6+N):(5+2*N)] , 2 , quantile, probs=0.05 )

Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[4001:40000,(6+N):(5+2*N)] , 2 , quantile, probs=0.95 )



Epars_m[1]=1
Epars_5[1]=1
Epars_9[1]=1

####################  trace plot



traceplot(fit, pars=c("u","b1","b3"), inc_warmup = T,  nrow = 3, ncol = 1, window = NULL, include = TRUE)



traceplot(fit, pars=c("lambda","sig2","sigob"), inc_warmup = T,  nrow = 3, ncol = 1, window = NULL, include = TRUE)



##################   plot 30/45/60/   length(y)

#  Germany Italy Japan Spain UK US Australia



m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main="In-sample fit results of M5 model for Germany",pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main="In-sample fit results of M5 model for Germany",pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main="In-sample fit results of M5 model for Germany",lwd = 2)












#############################################################
n=length(y)

u=0.1
b1=-0.0000005
b3=0


n=length(y)

x=c()
x[1]=100
for (i in 1:(n-1) ){

x[i+1]=((x[i]^2)*exp(u+b1*x[i]))/(b3+x[i])
}
x=round(x)
x


plot(y,ylim=c(min(x),max(x)))
par(new=TRUE)
plot(x,col="red",ylim=c(min(x),max(x)))















###############################################################################
## u=0.14,b1=-0.00000009, b3=5
# u ~ normal(0,1);
# b1 ~ normal(0,0.0000001);
# b3 ~ normal(0,10);







### Germany
u=0.14
b1=-0.0000009
b2=0.4



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


