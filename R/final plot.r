###############################################################################
# "Infection Rate Models for COVID-19: 
#  Model Risk and Public Health News Sentiment Exposure Adjustments"
#   
#  Ioannis Chalkiadakis, Kevin Hongxuan Yan, Gareth W. Peters, Pavel V. Shevchenko
#
#  Kevin Hongxuan Yan
#  March 2021
###############################################################################


####### read new data Y
setwd("./covid19modelrisk/data/covid19_infected/")

  y_GM=read.csv("COVID-19_confirmed__Germany.csv")[,2]
  y_IT=read.csv("COVID-19_confirmed__Italy.csv")[,2]
  y_JP=read.csv("COVID-19_confirmed__Japan.csv")[,2]
  y_SP=read.csv("COVID-19_confirmed__Spain.csv")[,2]
  y_UK=read.csv("COVID-19_confirmed__United Kingdom.csv")[,2]
  y_US=read.csv("COVID-19_confirmed__US.csv")[,2]
  y_AU=read.csv("COVID-19_confirmed_Australia.csv")[,2]



name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")

pdf("UK_tplot.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_UK,main=name[5],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,length(y_UK)), labels=c("Jan/31","May/09","Aug/06"))   ### UK
dev.off()


pdf("GM_tplot.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_GM,main=name[1],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,length(y_GM)), labels=c("Jan/27","May/05","Aug/06"))   ### GM
dev.off()


pdf("IT_tplot.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_IT,main=name[2],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,length(y_IT)), labels=c("Jan/31","May/09","Aug/06"))   ### IT
dev.off()


pdf("SP_tplot.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_SP,main=name[4],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,length(y_SP)), labels=c("Feb/01","May/10","Aug/06"))   ### SP
dev.off()


pdf("US_tplot.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_US,main=name[6],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,length(y_US)), labels=c("Jan/22","Apr/30","Aug/06"))   ### US
dev.off()


pdf("JP_tplot.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_JP,main=name[3],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,length(y_JP)), labels=c("Jan/22","Apr/30","Aug/06"))   ### JP
dev.off()


pdf("AU_tplot.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_AU,main=name[7],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,length(y_AU)), labels=c("Jan/26","May/04","Aug/06"))   ### AU
dev.off()






















####################  Violin plot 

library(methods)
library(rstan)


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")



k=7


setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\M2 Rdata")
load(paste("M2_",abn[k],".Rdata", sep=""))
Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , mean )
res=Epars_m-y
Res_2=(res-mean(res))/(sd(res))
R_2=cbind(Res_2,rep("M2",length(Res_2)))


setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\M4 Rdata")
load(paste("M4_",abn[k],".Rdata", sep=""))
Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , mean )
res=Epars_m-y
Res_4=(res-mean(res))/(sd(res))
R_4=cbind(Res_4,rep("M4",length(Res_4)))


setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\M7 Rdata")
load(paste("M7_",abn[k],".Rdata", sep=""))
Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , mean )
res=Epars_m-y
Res_7=(res-mean(res))/(sd(res))
R_7=cbind(Res_7,rep("M7",length(Res_7)))


setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\M8 Rdata")
load(paste("M8_",abn[k],".Rdata", sep=""))
Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , mean )
res=Epars_m-y
Res_8=(res-mean(res))/(sd(res))
R_8=cbind(Res_8,rep("M8",length(Res_8)))


setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\M9 Rdata")
load(paste("M9_",abn[k],".Rdata", sep=""))
Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , mean )
res=Epars_m-y
Res_9=(res-mean(res))/(sd(res))
R_9=cbind(Res_9,rep("M9",length(Res_9)))


setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\M10 Rdata")
load(paste("M10_",abn[k],".Rdata", sep=""))
Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , mean )
res=Epars_m-y
Res_10=(res-mean(res))/(sd(res))
R_10=cbind(Res_10,rep("M10",length(Res_10)))


setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\M11 Rdata")
load(paste("M11_",abn[k],".Rdata", sep=""))
Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(8+N):(7+2*N)] , 2 , mean )
res=Epars_m-y
Res_11=(res-mean(res))/(sd(res))
R_11=cbind(Res_11,rep("M11",length(Res_11)))


setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\M12 Rdata")
load(paste("M12_",abn[k],".Rdata", sep=""))
Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , mean )
res=Epars_m-y
Res_12=(res-mean(res))/(sd(res))
R_12=cbind(Res_12,rep("M12",length(Res_12)))



setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\final plot")

pdf(paste(abn[k],"_plot.pdf",sep=""), width = 16.5, height = 8.50) 

library("vioplot")
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
vioplot(Res_2,Res_4,Res_7,Res_8,Res_9,Res_10,Res_11,Res_12,names=c("M2","M4","M7","M8","M9","M10","M11","M12"))
title(paste("Residual plots for ",name[k],sep=""))


dev.off()



##############   boxplot   ########################################################################################################################################






library(methods)
library(rstan)


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")



k=3

setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\M2 Rdata")    ##### change
load(paste("M2_",abn[k],".Rdata", sep=""))    ##### change
Epars <- as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)]


res_all=matrix(,nrow=dim(Epars)[1],ncol=dim(Epars)[2])
Y_m=matrix(rep(y,nrow=dim(Epars)[1]),nrow=dim(Epars)[1],ncol=dim(Epars)[2],byrow=T)
res_all=as.matrix(Epars-Y_m)


res_st=matrix(,nrow=dim(Epars)[1],ncol=dim(Epars)[2])
for (i in 1:dim(Epars)[1]){
res_st[i,]=(res_all[i,]-mean(res_all[i,]))/(sd(res_all[i,]))
}





setwd("C:\\Users\\yhx19\\Desktop\\NEW RUN\\final plot")


pdf(paste(abn[k],"_boxplot.pdf",sep=""), width = 16.5, height = 8.50) 

par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
boxplot(res_st,xlab="time steps")
title(paste("M2 residual boxplots for ",name[k],sep=""))   ##### change


dev.off()








##############   NEW plot   ########################################################################################################################################
##############   NEW plot   ########################################################################################################################################



####### read new data Y
setwd("G:\\Covid paper\\additionalcountdataofinfectedcases\\new data")

  y_GM=read.csv("COVID-19_confirmed__Germany_latest_Feb.csv")[,2]
  y_IT=read.csv("COVID-19_confirmed__Italy_latest_Feb.csv")[,2]
  y_JP=read.csv("COVID-19_confirmed__Japan_latest_Feb.csv")[,2]
  y_SP=read.csv("COVID-19_confirmed__Spain_latest_Feb.csv")[,2]
  y_UK=read.csv("COVID-19_confirmed__United Kingdom_latest_Feb.csv")[,2]
  y_US=read.csv("COVID-19_confirmed__US_latest_Feb.csv")[,2]
  y_AU=read.csv("COVID-19_confirmed_Australia_latest_Feb.csv")[,2]


D_GM=c("Jan/27","May/05","Aug/13","Nov/21","Jan/10")
D_IT=c("Jan/31","May/09","Aug/17","Nov/25","Jan/14")
D_JP=c("Jan/22","Apr/30","Aug/08","Nov/16","Jan/5")
D_SP=c("Feb/01","May/10","Aug/18","Nov/26","Jan/15")
D_UK=c("Jan/31","May/09","Aug/17","Nov/25","Jan/14")
D_US=c("Jan/22","Apr/30","Aug/08","Nov/16","Jan/5")
D_AU=c("Jan/26","May/04","Aug/12","Nov/20","Jan/9")

Date=c(D_GM,D_IT,D_JP,D_SP,D_UK,D_US,D_AU)


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")

k=5
pdf("UK_tplot_LONG.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_UK,main=name[5],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,200,300,350), labels=c("Jan/31","May/09","Aug/17","Nov/25","Jan/14"))   
   ### UK
dev.off()


k=1
pdf("GM_tplot_LONG.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_GM,main=name[1],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,200,300,350), labels=c("Jan/27","May/05","Aug/13","Nov/21","Jan/10"))   
   ### GM
dev.off()

k=2
pdf("IT_tplot_LONG.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_IT,main=name[2],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,200,300,350), labels=c("Jan/31","May/09","Aug/17","Nov/25","Jan/14"))   
   ### IT
dev.off()

k=4
pdf("SP_tplot_LONG.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_SP,main=name[4],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,200,300,350), labels=c("Feb/01","May/10","Aug/18","Nov/26","Jan/15"))   
   ### SP
dev.off()

k=6
pdf("US_tplot_LONG.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_US,main=name[6],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,200,300,350), labels=c("Jan/22","Apr/30","Aug/08","Nov/16","Jan/5"))   
   ### US
dev.off()

k=3
pdf("JP_tplot_LONG.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_JP,main=name[3],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,200,300,350), labels=c("Jan/22","Apr/30","Aug/08","Nov/16","Jan/5"))   
   ### JP
dev.off()

k=7
pdf("AU_tplot_LONG.pdf", width = 16.5, height = 8.50) 
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y_AU,main=name[7],xlab="time",ylab="count",xaxt = 'n')
axis(1, at=c(0,100,200,300,350), labels=c("Jan/26","May/04","Aug/12","Nov/20","Jan/9"))   
   ### AU
dev.off()























