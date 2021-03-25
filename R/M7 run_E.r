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
u ~ normal(0.2,0.1);
b5 ~ normal(1200000,100000);
b6 ~ normal(1200000,100000);
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
  y=read.csv("COVID-19_confirmed__US.csv")
##  y=read.csv("COVID-19_confirmed_Australia.csv")

y=y[,2]

setwd("./covid19modelrisk/output/model_with_E")


###########
##  E=read.csv("volumeWeight_entropy_sentiment_IQRSummary_GM.csv")
##  E=read.csv("volumeWeight_entropy_sentiment_IQRSummary_IT.csv")
##  E=read.csv("volumeWeight_entropy_sentiment_IQRSummary_JP.csv")
##  E=read.csv("volumeWeight_entropy_sentiment_IQRSummary_SP.csv")
##  E=read.csv("volumeWeight_entropy_sentiment_IQRSummary_UK.csv")
  E=read.csv("volumeWeight_entropy_sentiment_IQRSummary_US.csv")
##  E=read.csv("volumeWeight_entropy_sentiment_IQRSummary_AU.csv")





################
  E=E[,4]
##  E=E[,5]



N=length(y)
dat <- list(N = N, y = y,E=E)


#### initial values
inits =list(list(lambda=0,u=0.2,b5=1200000, b6=1200000, sig2=0.01,e=rep(0.01,N),sigob=1000))  ### 134


fit <- stan(model_code = stanmodelcode, model_name = "example",init =inits,
data = dat, warmup=10000,iter = 100000, chains = 1, seed=999,thin=1)
print(fit,digits = 4)








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








