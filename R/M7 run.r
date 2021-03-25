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
real b5;
real b6;
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
	
      mu[n] <- mu[n-1]*exp(u*exp(-((mu[n-1]-b5)/b6)^2)+e[n]); 
	  
	  }

}



##### define piror and likelihood function
model {
lambda~uniform(-1,1);
u ~ normal(0.07,0.1);
b5 ~ normal(3000,100000);
b6 ~ normal(3000,100000);
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
setwd("./covid19modelrisk/output/forecast")

##  y=read.csv("COVID-19_confirmed__Germany.csv")
##  y=read.csv("COVID-19_confirmed__Italy.csv")
##  y=read.csv("COVID-19_confirmed__Japan.csv")
##  y=read.csv("COVID-19_confirmed__Spain.csv")
##  y=read.csv("COVID-19_confirmed__United Kingdom.csv")
##  y=read.csv("COVID-19_confirmed__US.csv")
  y=read.csv("COVID-19_confirmed_Australia.csv")

y=y[,2]




E=rep(1,length(y))

N=length(y)
dat <- list(N = N, y = y,E=E)


#### initial values
inits =list(list(lambda=0,u=0.07,b5=3000, b6=3000, sig2=0.01,e=rep(0.01,N),sigob=1000))  ### 134


fit <- stan(model_code = stanmodelcode, model_name = "example",init =inits,
data = dat, warmup=10000,iter = 100000, chains = 1, seed=999,thin=1)
print(fit,digits = 4)





####################  trace plot

library(methods)
library(rstan)


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")




k=6


pdf(paste("M7_",abn[k],"_p1.pdf",sep=""), width = 16.5, height = 8.50) 

traceplot(fit, pars=c("u","b5","b6"), inc_warmup = T,  nrow = 3, ncol = 1, window = NULL, include = TRUE)

dev.off()

pdf(paste("M7_",abn[k],"_p2.pdf",sep=""), width = 16.5, height = 8.50) 

traceplot(fit, pars=c("sigob","sig2","lambda"), inc_warmup = T,  nrow = 3, ncol = 1, window = NULL, include = TRUE)
dev.off()



##################   plot 30/45/60/   length(y)


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , mean )

Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.05 )

Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.95 )






pdf(paste("M7_",abn[k],"_30.pdf",sep=""), width = 16.5, height = 8.50) 

m=30

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M7_",abn[k],"_45.pdf",sep=""), width = 16.5, height = 8.50) 

m=45

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M7_",abn[k],"_60.pdf",sep=""), width = 16.5, height = 8.50) 

m=60

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),lwd = 2)

dev.off()





pdf(paste("M7_",abn[k],"_all.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),lwd = 2)

dev.off()






########    log plot          ####################################################

k=1


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")




Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , mean )

Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.05 )

Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.95 )




y=log(y)
Epars_m=log(Epars_m)
Epars_9=log(Epars_9)
Epars_5=log(Epars_5)


setwd("./covid19modelrisk/output/plot")

pdf(paste("M7_",abn[k],"_log.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19)
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),lwd = 2)

dev.off()





#####################################################################################################################################################
######################################      forecast          ###################################

################################################



P=20
Mu=matrix(,nrow=90000,ncol=P)
X=matrix(,nrow=90000,ncol=P)


lambda_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,1]
u_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,2]
b5_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,3]
b6_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,4]
sig2_p=as.data.frame(fit@sim[[1]][[1]])[10001:100000,5]


set.seed(999)

for (r in 1:90000){

Mu[r,1]=y[length(y)]
for (s in 1:(P-1) ){

e=rnorm(1,0,sig2_p[r])

Mu[r,s+1]=Mu[r,s]*exp(u_p[r]*exp(-((Mu[r,s]-b5_p[r])/b6_p[r])^2)+e)

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
setwd("./covid19modelrisk/output/forecast")

##  FT=read.csv("COVID-19_confirmed__Germany_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed__Italy_FT.csv")[,2]
FT=read.csv("COVID-19_confirmed__Japan_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed__Spain_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed__United Kingdom_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed__US_FT.csv")[,2]
##  FT=read.csv("COVID-19_confirmed_Australia_FT.csv")[,2]





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


#####
#######
############

k=3



name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")


m=length(Zm)


pdf(paste("M7_",abn[k],"_F.pdf",sep=""), width = 16.5, height = 8.50) 

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(Zm[1:m],ylim=c(min(Zm[1:m]),max(Zm[1:m],Z9[1:m],Z_FT)),ylab="Count",xlab="Time",main=paste("Out-sample forecast results of M7 model for ",name[k],sep=""),pch=19,col="red",xaxt = 'n')
polygon(c(rev(newx), newx), c(rev(Z9[1:m]), Z5[1:m]), col = "grey80", border = NA)

par(new=TRUE)
plot(Zm[1:m],ylim=c(min(Zm[1:m]),max(Zm[1:m],Z9[1:m],Z_FT)),ylab="Count",xlab="Time",main=paste("Out-sample forecast results of M7 model for ",name[k],sep=""),pch=19,col="red",xaxt = 'n')
lines(newx,Z9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Z5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(Z_FT[1:m],ylim=c(min(Zm[1:m]),max(Zm[1:m],Z9[1:m],Z_FT)),ylab="Count",xlab="Time",main=paste("Out-sample forecast results of M7 model for ",name[k],sep=""),pch=19,xaxt = 'n')
legend("topleft", legend=c("Observed data", "Forecasted data"), col=c("black", "red"), lty=c(1,1), lwd = 3,cex=2)

# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/27","Mar/17","May/05","Jun/25","Aug/06"))   ### GM
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/31","Mar/13","May/09","Jun/29","Aug/06"))   ### IT
axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/22","Mar/12","Apr/30","Jun/19","Aug/06"))   ### JP
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Feb/01","Mar/20","May/10","Jun/30","Aug/06"))   ### SP
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/31","Mar/21","May/09","Jun/29","Aug/06"))   ### UK
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/22","Mar/12","Apr/30","Jun/19","Aug/06"))   ### US
# axis(1, at=c(0,50,100,150,length(y)), labels=c("Jan/26","Mar/16","May/04","Jun/24","Aug/06"))   ### AU

dev.off()







######   RMSE results        ###############################################################

### Germany
> list(mae,rmse,mape,mase)
[[1]]
[1] 8027.949

[[2]]
[1] 9396.267

[[3]]
[1] 0.03476308

[[4]]
[1] 6.685852


## Italy
> list(mae,rmse,mape,mase)
[[1]]
[1] 7980.055

[[2]]
[1] 9472.923

[[3]]
[1] 0.03104283

[[4]]
[1] 11.86022


## Japan
> list(mae,rmse,mape,mase)
[[1]]
[1] 11774.83

[[2]]
[1] 13085.23

[[3]]
[1] 0.2015648

[[4]]
[1] 11.76554



## Spain
> list(mae,rmse,mape,mase)
[[1]]
[1] 41519.09

[[2]]
[1] 50731.13

[[3]]
[1] 0.1096009

[[4]]
[1] 7.478292



## UK
> list(mae,rmse,mape,mase)
[[1]]
[1] 8588.052

[[2]]
[1] 9873.951

[[3]]
[1] 0.02666235

[[4]]
[1] 8.224031



## US
> list(mae,rmse,mape,mase)
[[1]]
[1] 65821.77

[[2]]
[1] 77369.36

[[3]]
[1] 0.01208277

[[4]]
[1] 1.418199


## AU
> list(mae,rmse,mape,mase)
[[1]]
[1] 2985.154

[[2]]
[1] 3266.54

[[3]]
[1] 0.1244203

[[4]]
[1] 12.26599










###########################################################################################################################
###########################################################################################################################
###########################################################################################################################


mean(y)
median(y)

n=length(y)



u=0.07
b5=3000
b6=3000


x=c()
x[1]=100
for (i in 1:(n-1) ){

x[i+1]=x[i]*exp(u*exp(-((b5-x[i])/b6)^2))
}
x=round(x)
x



plot(y,ylim=c(min(y,x),max(y,x)))
par(new=TRUE)
plot(x,col="red",ylim=c(min(y,x),max(y,x)))






#############################################################
##### GM

u=0.25
b5=60000
b6=60000



##### IT

u=0.25
b5=70000
b6=70000


##### JP
u=0.1
b5=10000
b6=10000



##### SP
u=0.25
b5=70000
b6=70000



##### UK
u=0.2
b5=80000
b6=80000




##### US
u=0.2
b5=1200000
b6=1200000




##### AU
u=0.07
b5=3000
b6=3000






