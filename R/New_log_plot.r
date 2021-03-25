###############################################################################
# "Infection Rate Models for COVID-19: 
#  Model Risk and Public Health News Sentiment Exposure Adjustments"
#   
#  Ioannis Chalkiadakis, Kevin Hongxuan Yan, Gareth W. Peters, Pavel V. Shevchenko
#
#  Kevin Hongxuan Yan
#  March 2021
###############################################################################

for (k in 1:7){


D_GM=c("Jan/27","Mar/17","May/05","Jun/25","Jul/17")
D_IT=c("Jan/31","Mar/13","May/09","Jun/29","Jul/17")
D_JP=c("Jan/22","Mar/12","Apr/30","Jun/19","Jul/17")
D_SP=c("Feb/01","Mar/20","May/10","Jun/30","Jul/17")
D_UK=c("Jan/31","Mar/21","May/09","Jun/29","Jul/17")
D_US=c("Jan/22","Mar/12","Apr/30","Jun/19","Jul/17")
D_AU=c("Jan/26","Mar/16","May/04","Jun/24","Jul/17")

Date=c(D_GM,D_IT,D_JP,D_SP,D_UK,D_US,D_AU)


name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")



setwd("./covid19modelrisk/output/M4")

load(paste("M2_",abn[k],".RData",sep=""))


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , mean )
Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , quantile, probs=0.05 )
Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , quantile, probs=0.95 )


y=log(y)
Epars_m=log(Epars_m)
Epars_9=log(Epars_9)
Epars_5=log(Epars_5)

setwd("./covid19modelrisk/output/plot")

pdf(paste("M2_",abn[k],"_log.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),lwd = 2,xaxt = 'n')
legend("bottomright", legend=c("Observed data", "Fitted data"), col=c("black", "red"), lty=c(1,1), lwd = 3,cex=2)
axis(1, at=c(0,50,100,150,180), labels=c(Date[(k-1)*5+1],Date[(k-1)*5+2],Date[(k-1)*5+3],Date[(k-1)*5+4],Date[(k-1)*5+5]))   
dev.off()






setwd("./covid19modelrisk/output/M2")

load(paste("M4_",abn[k],".RData",sep=""))


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , mean )
Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.05 )
Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.95 )


y=log(y)
Epars_m=log(Epars_m)
Epars_9=log(Epars_9)
Epars_5=log(Epars_5)

setwd("./covid19modelrisk/output/plot")

pdf(paste("M4_",abn[k],"_log.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M4 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M4 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M4 model for ",name[k],sep=""),lwd = 2,xaxt = 'n')
legend("bottomright", legend=c("Observed data", "Fitted data"), col=c("black", "red"), lty=c(1,1), lwd = 3,cex=2)
axis(1, at=c(0,50,100,150,180), labels=c(Date[(k-1)*5+1],Date[(k-1)*5+2],Date[(k-1)*5+3],Date[(k-1)*5+4],Date[(k-1)*5+5]))   
dev.off()






setwd("./covid19modelrisk/output/M7")

load(paste("M7_",abn[k],".RData",sep=""))


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
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),lwd = 2,xaxt = 'n')
legend("bottomright", legend=c("Observed data", "Fitted data"), col=c("black", "red"), lty=c(1,1), lwd = 3,cex=2)
axis(1, at=c(0,50,100,150,180), labels=c(Date[(k-1)*5+1],Date[(k-1)*5+2],Date[(k-1)*5+3],Date[(k-1)*5+4],Date[(k-1)*5+5]))   
dev.off()





setwd("./covid19modelrisk/output/M12")

load(paste("M12_",abn[k],".RData",sep=""))


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , mean )
Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.05 )
Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.95 )


y=log(y)
Epars_m=log(Epars_m)
Epars_9=log(Epars_9)
Epars_5=log(Epars_5)

setwd("./covid19modelrisk/output/plot")

pdf(paste("M12_",abn[k],"_log.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M12 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M12 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M12 model for ",name[k],sep=""),lwd = 2,xaxt = 'n')
legend("bottomright", legend=c("Observed data", "Fitted data"), col=c("black", "red"), lty=c(1,1), lwd = 3,cex=2)
axis(1, at=c(0,50,100,150,180), labels=c(Date[(k-1)*5+1],Date[(k-1)*5+2],Date[(k-1)*5+3],Date[(k-1)*5+4],Date[(k-1)*5+5]))   
dev.off()





D_GM=c("Jan/27","May/05","Aug/13","Nov/21","Jan/10")
D_IT=c("Jan/31","May/09","Aug/17","Nov/25","Jan/14")
D_JP=c("Jan/22","Apr/30","Aug/08","Nov/16","Jan/5")
D_SP=c("Feb/01","May/10","Aug/18","Nov/26","Jan/15")
D_UK=c("Jan/31","May/09","Aug/17","Nov/25","Jan/14")
D_US=c("Jan/22","Apr/30","Aug/08","Nov/16","Jan/5")
D_AU=c("Jan/26","May/04","Aug/12","Nov/20","Jan/9")

Date=c(D_GM,D_IT,D_JP,D_SP,D_UK,D_US,D_AU)



name=c("Germany","Italy","Japan","Spain","U.K.","U.S.","Australia")
abn=c("GM","IT","JP","SP","UK","US","AU")



setwd("./covid19modelrisk/output/M2_extra")

load(paste("M2_",abn[k],".RData",sep=""))


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , mean )
Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , quantile, probs=0.05 )
Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(6+N):(5+2*N)] , 2 , quantile, probs=0.95 )


y=log(y)
Epars_m=log(Epars_m)
Epars_9=log(Epars_9)
Epars_5=log(Epars_5)

setwd("./covid19modelrisk/output/plot")

pdf(paste("M2_",abn[k],"_log_LONG.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M2 model for ",name[k],sep=""),lwd = 2,xaxt = 'n')
legend("bottomright", legend=c("Observed data", "Fitted data"), col=c("black", "red"), lty=c(1,1), lwd = 3,cex=2)
axis(1, at=c(0,100,200,300,350), labels=c(Date[(k-1)*5+1],Date[(k-1)*5+2],Date[(k-1)*5+3],Date[(k-1)*5+4],Date[(k-1)*5+5]))   
  
dev.off()









setwd("./covid19modelrisk/output/M7_extra")

load(paste("M7_",abn[k],".RData",sep=""))


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , mean )
Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.05 )
Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.95 )


y=log(y)
Epars_m=log(Epars_m)
Epars_9=log(Epars_9)
Epars_5=log(Epars_5)

setwd("./covid19modelrisk/output/plot")

pdf(paste("M7_",abn[k],"_log_LONG.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M7 model for ",name[k],sep=""),lwd = 2,xaxt = 'n')
legend("bottomright", legend=c("Observed data", "Fitted data"), col=c("black", "red"), lty=c(1,1), lwd = 3,cex=2)
axis(1, at=c(0,100,200,300,350), labels=c(Date[(k-1)*5+1],Date[(k-1)*5+2],Date[(k-1)*5+3],Date[(k-1)*5+4],Date[(k-1)*5+5]))   
   
dev.off()





setwd("./covid19modelrisk/output/M12_extra")

load(paste("M12_",abn[k],".RData",sep=""))


Epars_m <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , mean )
Epars_5 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.05 )
Epars_9 <- apply( as.data.frame(fit@sim[[1]][[1]])[10001:100000,(7+N):(6+2*N)] , 2 , quantile, probs=0.95 )


y=log(y)
Epars_m=log(Epars_m)
Epars_9=log(Epars_9)
Epars_5=log(Epars_5)

setwd("./covid19modelrisk/output/plot")

pdf(paste("M12_",abn[k],"_log_LONG.pdf",sep=""), width = 16.5, height = 8.50) 

m=length(y)

newx=c(1:m)
par(cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mar=c(5, 5, 3, 0.5))
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M12 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)

polygon(c(rev(newx), newx), c(rev(Epars_9[1:m]), Epars_5[1:m]), col = "grey80", border = NA)

lines(newx,Epars_9[1:m], lty = 'dashed', col = 'red',lwd = 2)
lines(newx, Epars_5[1:m], lty = 'dashed', col = 'red',lwd = 2)
par(new=TRUE)
plot(y[1:m],ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M12 model for ",name[k],sep=""),pch=19,xaxt = 'n')
par(new=TRUE)
plot.ts(Epars_m[1:m],col="red",ylim=c(min(y[1:m]),max(y[1:m],Epars_9[1:m])),ylab="Count",xlab="Time",main=paste("In-sample fit results of M12 model for ",name[k],sep=""),lwd = 2,xaxt = 'n')
legend("bottomright", legend=c("Observed data", "Fitted data"), col=c("black", "red"), lty=c(1,1), lwd = 3,cex=2)
axis(1, at=c(0,100,200,300,350), labels=c(Date[(k-1)*5+1],Date[(k-1)*5+2],Date[(k-1)*5+3],Date[(k-1)*5+4],Date[(k-1)*5+5]))   
   
dev.off()





}





















