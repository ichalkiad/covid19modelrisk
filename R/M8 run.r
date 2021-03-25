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
		  ## +log(1-normal_cdf(0.08*u,u,0.01))*(b0>0.08*u)
		  ## +log(normal_cdf((0.515+0.125*u),u,0.01))*(b2>(0.515+0.125*u));    #### +log(b0>0.08*u); (b0>0.08*u)
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
real b1;
real b2;
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
	
      mu[n] <- mu[n-1]*exp(u*exp(-b0*n)+b1*(mu[n-1]^b2)+e[n]); 
	  
	  }

}



##### define piror and likelihood function
model {
lambda~uniform(-1,1);
u ~ normal(0.13,0.1);
b0 ~ normal(0.025,0.01);
b1 ~ normal(-0.00000006,0.00000001);
b2 ~ normal(0.1,0.1);
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
setwd("./covid19modelrisk/data/covid19_infected")

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
inits =list(list(lambda=0,u=0.1,b0=0.025,b1=-0.00000006, b2=0.1,sig2=0.0005,e=rep(0.01,N),sigob=1000),
list(lambda=0,u=0.2,b0=0.025,b1=-0.00000006, b2=0.1,sig2=0.0005,e=rep(0.01,N),sigob=1000),
list(lambda=0,u=0.15,b0=0.025,b1=-0.00000006, b2=0.1,sig2=0.0005,e=rep(0.01,N),sigob=1000),
list(lambda=0,u=0.25,b0=0.025,b1=-0.00000006, b2=0.1,sig2=0.0005,e=rep(0.01,N),sigob=1000),
list(lambda=0,u=0.2,b0=0.05,b1=-0.00000006, b2=0.1,sig2=0.0005,e=rep(0.01,N),sigob=1000),
list(lambda=0,u=0.2,b0=0.01,b1=-0.00000006, b2=0.1,sig2=0.0005,e=rep(0.01,N),sigob=1000),
list(lambda=0,u=0.2,b0=0.01,b1=-0.0000001, b2=0.1,sig2=0.0005,e=rep(0.01,N),sigob=1000),
list(lambda=0,u=0.2,b0=0.01,b1=-0.00000001, b2=0.1,sig2=0.0005,e=rep(0.01,N),sigob=1000))


fit <- stan(model_code = stanmodelcode, model_name = "example",init =inits,
data = dat, warmup=10000,iter = 100000, chains = 8, seed=999,thin=1)
print(fit,digits = 4)










####################  trace plot

library(methods)
library(rstan)


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")




k=1


pdf(paste("M8_",abn[k],"_p1.pdf",sep=""), width = 16.5, height = 8.50) 

traceplot(fit, pars=c("u","b0","b1","b2"), inc_warmup = T,  nrow = 4, ncol = 1, window = NULL, include = TRUE)

dev.off()

pdf(paste("M8_",abn[k],"_p2.pdf",sep=""), width = 16.5, height = 8.50) 

traceplot(fit, pars=c("sigob","sig2","lambda"), inc_warmup = T,  nrow = 3, ncol = 1, window = NULL, include = TRUE)
dev.off()



##################   plot 30/45/60/   length(y)


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , mean )

Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , quantile, probs=0.05 )

Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , quantile, probs=0.95 )




pdf(paste("M8_",abn[k],"_30.pdf",sep=""), width = 16.5, height = 8.50) 

m=30

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M8_",abn[k],"_45.pdf",sep=""), width = 16.5, height = 8.50) 

m=45

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M8_",abn[k],"_60.pdf",sep=""), width = 16.5, height = 8.50) 

m=60

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M8_",abn[k],"_all.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),lwd = 2)

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

pdf(paste("M8_",abn[k],"_log.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M8 model for ",name[k],sep=""),lwd = 2)

dev.off()







###########################################################################################################################





mean(y)
median(y)

n=length(y)



u=0.13
b0=0.025
b1=-0.00000006
b2=0.1



x=c()
x[1]=100
for (i in 1:(n-1) ){

x[i+1]=x[i]*exp(u*exp(-b0*i)+b1*(x[i]^b2))
}
x=round(x)
x



plot(y,ylim=c(min(y,x),max(y,x)))
par(new=TRUE)
plot(x,col="red",ylim=c(min(y,x),max(y,x)))






#############################################################
##### GM


u=0.2
b0=0.025
b1=-0.00000006
b2=0.1





##### IT

u=0.2
b0=0.025
b1=-0.00000006
b2=0.1



##### JP


u=0.15
b0=0.025
b1=-0.00000006
b2=0.1



##### SP

u=0.2
b0=0.025
b1=-0.00000006
b2=0.1




##### UK

u=0.2
b0=0.025
b1=-0.00000006
b2=0.1





##### US


u=0.27
b0=0.025
b1=-0.00000006
b2=0.1



##### AU


u=0.13
b0=0.025
b1=-0.00000006
b2=0.1












