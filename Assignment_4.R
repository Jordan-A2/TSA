## Assignment 4

## Question 4.1
data <- data.frame(read.csv("A4_hcab.csv"))
train.data<-data[1:727,1:3]
test.data<-data[728:751,1:3]
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(train.data[,2],type="l",col="red",xlim=c(1,751),ylim=c(4.6,404.53),
     main="Concentrations of the two gasses",xlab="Observation",
     ylab="Concentration (micrograms/cubic meter)")
lines(train.data[,3],col="blue")
lines(test.data[,2]~c(728:751),col="red",lty=2)
lines(test.data[,3]~c(728:751),col="blue",lty=2)
legend("top",legend=c("NO2","NOx","Training data","Testing data"),
       col=c("red","blue","black","black"),lty=c(1,1,1,2), cex=0.8)

NO2.diff <- rep(NA,559)
for (i in 1:559){
  NO2.diff[i] <- train.data[i+168,2]-train.data[i,2]
}

NOx.diff <- rep(NA,559)
for (i in 1:559){
  NOx.diff[i] <- train.data[i+168,3]-train.data[i,3]
}

plot(NO2.diff,type="l",col="red",xlim=c(1,559),ylim=c(-222.06,236.365),
     main="Weekly differences of the gases",xlab="Observation",
     ylab="Concentration (micrograms/cubic meter)")
lines(NOx.diff,col="blue")
legend("topleft",legend=c("NO2","NOx"),col=c("red","blue"),lty=1:1,cex=0.8)

## Question 4.3
initstate <- c(0,0)
initvar <- matrix(c(40,0,0,40),nrow=2)
A <- matrix(c(1,0,0,1),nrow=2)
C <- matrix(c(1,0,1,1),nrow=2)
sigma1 <- matrix(c(200,0,0,200),nrow=2)
sigma2 <- matrix(c(25,0,0,25),nrow=2)
library("FKF")
a<-fkf(a0=initstate,P0=initvar,dt=c(0,0),ct=c(0,0),Tt=A,Zt=C,
       HHt=sigma1,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
plot(a$att[1,],type="l",col="red",xlim=c(1,559),ylim=c(-135.7467,95.7729),
     main="Estimations of the weekly differences of NO",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")

plot(a$att[2,],type="l",col="red",xlim=c(1,559),ylim=c(-198.152,207.5478),
     main="Estimations of the weekly differences of NO2",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")


confintNO<-numeric(559)
confintNO2<-numeric(559)
for (i in 1:559){
  confintNO[i]<-qt(0.975,597)*sqrt(a$Ptt[1,1,i])
  confintNO2[i]<-qt(0.975,597)*sqrt(a$Ptt[2,2,i])
}

confintNO.lower<-a$att[1,]-confintNO
confintNO.upper<-a$att[1,]+confintNO
confintNO2.lower<-a$att[2,]-confintNO2
confintNO2.upper<-a$att[2,]+confintNO2

plot(a$att[1,],type="l",col="red",xlim=c(1,559),ylim=c(-148.1862,108.2124),
     main="95% confidence interval for the estimates of NO",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(confintNO.lower,col="green",lty=2)
lines(confintNO.upper,col="green",lty=2)
legend("bottomright",legend=c("95% confidence interval","Estimations"),
       col=c("green","red"),lty=c(2,1),cex=0.8)

plot(a$att[2,],type="l",col="red",xlim=c(1,559),ylim=c(-207.145,216.5409),
     main="95% confidence interval for the estimates of NO2",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(confintNO2.lower,col="green",lty=2)
lines(confintNO2.upper,col="green",lty=2)
legend("bottomright",legend=c("95% confidence interval","Estimations"),
       col=c("green","red"),lty=c(2,1),cex=0.8)

a$logLik

fcdiffstates <- matrix(NA,2,25)
fcsigma1 <- rep(list(matrix(NA,2,2)),25)
for (i in 1:24){
  fcdiffstates[,1] <- t(t(a$at[,560]))
  fcdiffstates[,i+1] <- A%*%fcdiffstates[,i]
}
for (i in 1:24){
  fcsigma1[[1]] <- a$Pt[,,560]
  fcsigma1[[i+1]] <- A%*%fcsigma1[[i]]%*%t(A)+sigma1
}      

fcdiffstates<-fcdiffstates[,2:25]
fcsigma1 <- fcsigma1[2:25]
fcdiff<-C%*%fcdiffstates
fcsigma2 <- rep(list(matrix(NA,2,2)),25)
for (i in 1:24){
  fcsigma2[[1]] <- fcsigma1[[1]]
  fcsigma2[[i+1]] <- C%*%fcsigma1[[i]]%*%t(C)+sigma2
} 
fcsigma2 <- fcsigma2[2:25]

fcNOx <- numeric(24)
fcNO2 <- numeric(24)
predintdiffNOx.lower <- numeric(24)
predintdiffNOx.upper <- numeric(24) 
predintdiffNO2.lower <- numeric(24)
predintdiffNO2.upper <- numeric(24)

for (i in 1:24){
predintdiffNOx.lower[i] <- fcdiff[1,i]-qt(0.975,597)*
  sqrt(a$Pt[1,1,560]+fcsigma2[[i]][1,1])
predintdiffNOx.upper[i] <- fcdiff[1,i]+qt(0.975,597)*
  sqrt(a$Pt[1,1,560]+fcsigma2[[i]][1,1])
predintdiffNO2.lower[i] <- fcdiff[2,i]-qt(0.975,597)*
  sqrt(a$Pt[2,2,560]+fcsigma2[[i]][2,2])
predintdiffNO2.upper[i] <- fcdiff[2,i]+qt(0.975,597)*
  sqrt(a$Pt[2,2,560]+fcsigma2[[i]][2,2])
}

predintNOx.lower <- numeric(24)
predintNOx.upper <- numeric(24)
predintNO2.lower <- numeric(24)
predintNO2.upper <- numeric(24)

for (i in 560:583){
  fcNOx[i-559] <- fcdiff[1,i-559]+train.data[i,3]
  fcNO2[i-559] <- fcdiff[2,i-559]+train.data[i,2]
  predintNOx.lower[i-559] <- predintdiffNOx.lower[i-559]+train.data[i,3]
  predintNOx.upper[i-559] <- predintdiffNOx.upper[i-559]+train.data[i,3]
  predintNO2.lower[i-559] <- predintdiffNO2.lower[i-559]+train.data[i,2]
  predintNO2.upper[i-559] <- predintdiffNO2.upper[i-559]+train.data[i,2]
}

plot(test.data[,3],type="l",col="red",xlim=c(1,24),ylim=c(-191.4945,472.1747),
     main="Predictions for the final 24 hours for NOx",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(fcNOx,col="blue")
lines(predintNOx.lower,col="green",lty=2)
lines(predintNOx.upper,col="green",lty=2)
legend("topleft",legend=c("95% prediction interval","Predictions",
                          "Testing data"),
       col=c("green","blue","red"),lty=c(2,1,1),cex=0.8)

plot(test.data[,2],type="l",col="red",xlim=c(1,24),ylim=c(-116.371,243.3492),
     main="Predictions for the final 24 hours for NO2",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(fcNO2,col="blue")
lines(predintNO2.lower,col="green",lty=2)
lines(predintNO2.upper,col="green",lty=2)
legend("topleft",legend=c("95% prediction interval","Predictions",
                          "Testing data"),
       col=c("green","blue","red"),lty=c(2,1,1),cex=0.8)

predict.table.a <- data.frame(test.data[,3],fcNOx,predintNOx.lower,
                            predintNOx.upper)
names(predict.table.a) <- c("Actual NOx concentrations",
                         "NOx predicted concentrations",
                         "NOx prediction interval lower bound",
                         "NOx prediction interval upper bound")
predict.table.a <- predict.table.a[c(1,2,6,24),]

predict.table.b<-data.frame(test.data[,2],fcNO2,predintNO2.lower,
                            predintNO2.upper)
names(predict.table.b)<-c("Actual NO2 concentrations",
                          "NO2 predicted concentrations",
                          "NO2 prediction interval lower bound",
                          "NO2 prediction interval upper bound")
predict.table.b <- predict.table.b[c(1,2,6,24),]

## Question 4.4
func1 <- function(par){
  initvaropt<-par*matrix(c(1,0,0,1),nrow=2)
  a1<-fkf(a0=initstate,P0=initvaropt,dt=c(0,0),ct=c(0,0),Tt=A,Zt=C,
         HHt=sigma1,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a1$logLik))
}

opt1 <- optim(40,func1,method = "L-BFGS-B")
initvaropt<-opt1$par*matrix(c(1,0,0,1),nrow=2)

func2 <- function(par){
  sigma1opt<-par*matrix(c(1,0,0,1),nrow=2)
  a2<-fkf(a0=initstate,P0=initvaropt,dt=c(0,0),ct=c(0,0),Tt=A,Zt=C,
          HHt=sigma1opt,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a2$logLik))
}

opt2 <- optim(200,func2,method = "L-BFGS-B")
sigma1opt<-opt2$par*matrix(c(1,0,0,1),nrow=2)

func3 <- function(par){
  sigma2opt<-par*matrix(c(1,0,0,1),nrow=2)
  a3<-fkf(a0=initstate,P0=initvaropt,dt=c(0,0),ct=c(0,0),Tt=A,Zt=C,
          HHt=sigma1opt,GGt=sigma2opt,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a3$logLik))
}

opt3 <- optim(25,func3,method = "L-BFGS-B")
sigma2opt<-opt3$par*matrix(c(1,0,0,1),nrow=2)



aopt<-fkf(a0=initstate,P0=initvaropt,dt=c(0,0),ct=c(0,0),Tt=A,Zt=C,
       HHt=sigma1opt,GGt=sigma2opt,
       yt=matrix(c(NOx.diff,NO2.diff),nrow=2))

aopt$logLik

fcdiffstatesopt <- matrix(NA,2,25)
fcsigma1opt <- rep(list(matrix(NA,2,2)),25)
for (i in 1:24){
  fcdiffstatesopt[,1] <- t(t(aopt$at[,560]))
  fcdiffstatesopt[,i+1] <- A%*%fcdiffstatesopt[,i]
}
for (i in 1:24){
  fcsigma1opt[[1]] <- aopt$Pt[,,560]
  fcsigma1opt[[i+1]] <- A%*%fcsigma1opt[[i]]%*%t(A)+sigma1opt
}      

fcdiffstatesopt<-fcdiffstatesopt[,2:25]
fcsigma1opt <- fcsigma1opt[2:25]
fcdiffopt<-C%*%fcdiffstatesopt
fcsigma2opt <- rep(list(matrix(NA,2,2)),25)
for (i in 1:24){
  fcsigma2opt[[1]] <- fcsigma1opt[[1]]
  fcsigma2opt[[i+1]] <- C%*%fcsigma1opt[[i]]%*%t(C)+sigma2opt
} 
fcsigma2opt <- fcsigma2opt[2:25]

fcNOxopt <- numeric(24)
fcNO2opt <- numeric(24)
predintdiffNOxopt.lower <- numeric(24)
predintdiffNOxopt.upper <- numeric(24) 
predintdiffNO2opt.lower <- numeric(24)
predintdiffNO2opt.upper <- numeric(24)

for (i in 1:24){
  predintdiffNOxopt.lower[i] <- fcdiffopt[1,i]-qt(0.975,597)*
    sqrt(aopt$Pt[1,1,560]+fcsigma2opt[[i]][1,1])
  predintdiffNOxopt.upper[i] <- fcdiffopt[1,i]+qt(0.975,597)*
    sqrt(aopt$Pt[1,1,560]+fcsigma2opt[[i]][1,1])
  predintdiffNO2opt.lower[i] <- fcdiffopt[2,i]-qt(0.975,597)*
    sqrt(aopt$Pt[2,2,560]+fcsigma2opt[[i]][2,2])
  predintdiffNO2opt.upper[i] <- fcdiffopt[2,i]+qt(0.975,597)*
    sqrt(aopt$Pt[2,2,560]+fcsigma2opt[[i]][2,2])
}

predintNOxopt.lower <- numeric(24)
predintNOxopt.upper <- numeric(24)
predintNO2opt.lower <- numeric(24)
predintNO2opt.upper <- numeric(24)

for (i in 560:583){
  fcNOxopt[i-559] <- fcdiffopt[1,i-559]+train.data[i,3]
  fcNO2opt[i-559] <- fcdiffopt[2,i-559]+train.data[i,2]
  predintNOxopt.lower[i-559] <- predintdiffNOxopt.lower[i-559]+train.data[i,3]
  predintNOxopt.upper[i-559] <- predintdiffNOxopt.upper[i-559]+train.data[i,3]
  predintNO2opt.lower[i-559] <- predintdiffNO2opt.lower[i-559]+train.data[i,2]
  predintNO2opt.upper[i-559] <- predintdiffNO2opt.upper[i-559]+train.data[i,2]
}

plot(test.data[,3],type="l",col="red",xlim=c(1,24),ylim=c(-443.771,746.0131),
     main="Predictions for the final 24 hours for NOx",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(fcNOxopt,col="blue")
lines(predintNOxopt.lower,col="green",lty=2)
lines(predintNOxopt.upper,col="green",lty=2)
legend("topleft",legend=c("95% prediction interval for optimised predictions",
                          "Optimised predictions","Testing data"),
       col=c("green","blue","red"),lty=c(2,1,1),cex=0.8)

plot(test.data[,2],type="l",col="red",xlim=c(1,24),ylim=c(-297.073,439.1049),
     main="Predictions for the final 24 hours for NO2",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(fcNO2opt,col="blue")
lines(predintNO2opt.lower,col="green",lty=2)
lines(predintNO2opt.upper,col="green",lty=2)
legend("topleft",legend=c("95% prediction interval for optimised predictions",
                          "Optimised predictions","Testing data"),
       col=c("green","blue","red"),lty=c(2,1,1),cex=0.8)

predict.table.opt.a <- data.frame(test.data[,3],fcNOxopt,predintNOxopt.lower,
                              predintNOxopt.upper)
names(predict.table.opt.a) <- c("Actual NOx concentrations",
                            "NOx predicted concentrations",
                            "NOx prediction interval lower bound",
                            "NOx prediction interval upper bound")
predict.table.opt.a <- predict.table.opt.a[c(1,2,6,24),]

predict.table.opt.b<-data.frame(test.data[,2],fcNO2opt,predintNO2opt.lower,
                            predintNO2opt.upper)
names(predict.table.opt.b)<-c("Actual NO2 concentrations",
                          "NO2 predicted concentrations",
                          "NO2 prediction interval lower bound",
                          "NO2 prediction interval upper bound")
predict.table.opt.b <- predict.table.opt.b[c(1,2,6,24),]


## Question 4.5
func5.1 <- function(par){
  Aopt1<-par*matrix(c(1,0),nrow=2)
  Aopt2<-par*matrix(c(0,1),nrow=2)
  Aopt<-cbind(Aopt1,Aopt2)
  a5.1<-fkf(a0=initstate,P0=initvar,dt=c(0,0),ct=c(0,0),Tt=Aopt,Zt=C,
          HHt=sigma1,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a5.1$logLik))
}

opt5.1 <- optim(c(1,1),func5.1,method = "L-BFGS-B")
Aopt<-matrix(c(opt5.1$par[1],0,0,opt5.1$par[2]),nrow=2)

func5.2 <- function(par){
  initvaropt2<-par*matrix(c(1,0,0,1),nrow=2)
  a5.2<-fkf(a0=initstate,P0=initvaropt2,dt=c(0,0),ct=c(0,0),Tt=Aopt,Zt=C,
          HHt=sigma1,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a5.2$logLik))
}

opt5.2 <- optim(40,func5.2,method = "L-BFGS-B")
initvaropt2<-opt5.2$par*matrix(c(1,0,0,1),nrow=2)

func5.3 <- function(par){
  sigma1opt2<-par*matrix(c(1,0,0,1),nrow=2)
  a5.3<-fkf(a0=initstate,P0=initvaropt2,dt=c(0,0),ct=c(0,0),Tt=Aopt,Zt=C,
          HHt=sigma1opt2,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a5.3$logLik))
}

opt5.3 <- optim(200,func5.3,method = "L-BFGS-B")
sigma1opt2<-opt5.3$par*matrix(c(1,0,0,1),nrow=2)

func5.4 <- function(par){
  sigma2opt2<-par*matrix(c(1,0,0,1),nrow=2)
  a5.4<-fkf(a0=initstate,P0=initvaropt2,dt=c(0,0),ct=c(0,0),Tt=Aopt,Zt=C,
          HHt=sigma1opt2,GGt=sigma2opt2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a5.4$logLik))
}

opt5.4 <- optim(25,func5.4,method = "L-BFGS-B")
sigma2opt2<-opt5.4$par*matrix(c(1,0,0,1),nrow=2)



aopt2<-fkf(a0=initstate,P0=initvaropt2,dt=c(0,0),ct=c(0,0),Tt=Aopt,Zt=C,
          HHt=sigma1opt2,GGt=sigma2opt2,
          yt=matrix(c(NOx.diff,NO2.diff),nrow=2))

aopt2$logLik

fcdiffstatesopt2 <- matrix(NA,2,25)
fcsigma1opt2 <- rep(list(matrix(NA,2,2)),25)
for (i in 1:24){
  fcdiffstatesopt2[,1] <- t(t(aopt2$at[,560]))
  fcdiffstatesopt2[,i+1] <- Aopt%*%fcdiffstatesopt2[,i]
}
for (i in 1:24){
  fcsigma1opt2[[1]] <- aopt2$Pt[,,560]
  fcsigma1opt2[[i+1]] <- Aopt%*%fcsigma1opt2[[i]]%*%t(Aopt)+sigma1opt2
}      

fcdiffstatesopt2<-fcdiffstatesopt2[,2:25]
fcsigma1opt2 <- fcsigma1opt2[2:25]
fcdiffopt2<-C%*%fcdiffstatesopt2
fcsigma2opt2 <- rep(list(matrix(NA,2,2)),25)
for (i in 1:24){
  fcsigma2opt2[[1]] <- fcsigma1opt2[[1]]
  fcsigma2opt2[[i+1]] <- C%*%fcsigma1opt2[[i]]%*%t(C)+sigma2opt2
} 
fcsigma2opt2 <- fcsigma2opt2[2:25]

fcNOxopt2 <- numeric(24)
fcNO2opt2 <- numeric(24)
predintdiffNOxopt2.lower <- numeric(24)
predintdiffNOxopt2.upper <- numeric(24) 
predintdiffNO2opt2.lower <- numeric(24)
predintdiffNO2opt2.upper <- numeric(24)

for (i in 1:24){
  predintdiffNOxopt2.lower[i] <- fcdiffopt2[1,i]-qt(0.975,597)*
    sqrt(aopt2$Pt[1,1,560]+fcsigma2opt2[[i]][1,1])
  predintdiffNOxopt2.upper[i] <- fcdiffopt2[1,i]+qt(0.975,597)*
    sqrt(aopt2$Pt[1,1,560]+fcsigma2opt2[[i]][1,1])
  predintdiffNO2opt2.lower[i] <- fcdiffopt2[2,i]-qt(0.975,597)*
    sqrt(aopt2$Pt[2,2,560]+fcsigma2opt2[[i]][2,2])
  predintdiffNO2opt2.upper[i] <- fcdiffopt2[2,i]+qt(0.975,597)*
    sqrt(aopt2$Pt[2,2,560]+fcsigma2opt2[[i]][2,2])
}

predintNOxopt2.lower <- numeric(24)
predintNOxopt2.upper <- numeric(24)
predintNO2opt2.lower <- numeric(24)
predintNO2opt2.upper <- numeric(24)

for (i in 560:583){
  fcNOxopt2[i-559] <- fcdiffopt2[1,i-559]+train.data[i,3]
  fcNO2opt2[i-559] <- fcdiffopt2[2,i-559]+train.data[i,2]
  predintNOxopt2.lower[i-559] <- predintdiffNOxopt2.lower[i-559]+train.data[i,3]
  predintNOxopt2.upper[i-559] <- predintdiffNOxopt2.upper[i-559]+train.data[i,3]
  predintNO2opt2.lower[i-559] <- predintdiffNO2opt2.lower[i-559]+train.data[i,2]
  predintNO2opt2.upper[i-559] <- predintdiffNO2opt2.upper[i-559]+train.data[i,2]
}

plot(test.data[,3],type="l",col="red",xlim=c(1,24),ylim=c(-93.8527,421.6787),
     main="Predictions for the final 24 hours for NOx",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(fcNOxopt2,col="blue")
lines(predintNOxopt2.lower,col="green",lty=2)
lines(predintNOxopt2.upper,col="green",lty=2)
legend("topleft",legend=c("95% prediction interval for optimised predictions",
                          "Optimised predictions","Testing data"),
       col=c("green","blue","red"),lty=c(2,1,1),cex=0.8)

plot(test.data[,2],type="l",col="red",xlim=c(1,24),ylim=c(-84.25724,250),
     main="Predictions for the final 24 hours for NO2",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(fcNO2opt2,col="blue")
lines(predintNO2opt2.lower,col="green",lty=2)
lines(predintNO2opt2.upper,col="green",lty=2)
legend("topleft",legend=c("95% prediction interval for optimised predictions",
                          "Optimised predictions","Testing data"),
       col=c("green","blue","red"),lty=c(2,1,1),cex=0.8)

predict.table.opt2.a <- data.frame(test.data[,3],fcNOxopt2,predintNOxopt2.lower,
                                  predintNOxopt2.upper)
names(predict.table.opt2.a) <- c("Actual NOx concentrations",
                                "NOx predicted concentrations",
                                "NOx prediction interval lower bound",
                                "NOx prediction interval upper bound")
predict.table.opt2.a <- predict.table.opt2.a[c(1,2,6,24),]

predict.table.opt2.b<-data.frame(test.data[,2],fcNO2opt2,predintNO2opt2.lower,
                                predintNO2opt2.upper)
names(predict.table.opt2.b)<-c("Actual NO2 concentrations",
                              "NO2 predicted concentrations",
                              "NO2 prediction interval lower bound",
                              "NO2 prediction interval upper bound")
predict.table.opt2.b <- predict.table.opt2.b[c(1,2,6,24),]


## Question 4.6
func6.1a <- function(par){
  sigma1opt3.1<-par*matrix(c(1,0,0,0),nrow=2)
  sigma1opt3.2<-par*matrix(c(0,0,0,1),nrow=2)
  sigma1opt3a<-sigma1opt3.1+sigma1opt3.2
  a6.1a<-fkf(a0=initstate,P0=initvar,dt=c(0,0),ct=c(0,0),Tt=A,Zt=C,
            HHt=sigma1opt3a,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a6.1a$logLik))
}

opt6.1a <- optim(c(200,200),func6.1a,method = "L-BFGS-B")

func6.1b <- function(par){
  bbb <- par*c(1,1)
  sigma1opt3b<-matrix(c(1336.8185,bbb,987.0612),nrow=2)
  a6.1b<-fkf(a0=initstate,P0=initvar,dt=c(0,0),ct=c(0,0),Tt=A,Zt=C,
            HHt=sigma1opt3b,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a6.1b$logLik))
}

opt6.1b <- optim(0,func6.1b,method = "L-BFGS-B")
sigma1opt3<-matrix(c(opt6.1a$par[1],opt6.1b$par,
                     opt6.1b$par,opt6.1a$par[2]),nrow=2)

func6.2 <- function(par){
  Aopt4<-par*matrix(c(1,0),nrow=2)
  Aopt5<-par*matrix(c(0,1),nrow=2)
  Aopt3<-cbind(Aopt4,Aopt5)
  a6.2<-fkf(a0=initstate,P0=initvar,dt=c(0,0),ct=c(0,0),Tt=Aopt3,Zt=C,
            HHt=sigma1opt3,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a6.2$logLik))
}

opt6.2 <- optim(c(1,1),func6.2,method = "L-BFGS-B")
Aopt3<-matrix(c(opt6.2$par[1],0,0,opt6.2$par[2]),nrow=2)

func6.3 <- function(par){
  initvaropt3<-par*matrix(c(1,0,0,1),nrow=2)
  a6.3<-fkf(a0=initstate,P0=initvaropt3,dt=c(0,0),ct=c(0,0),Tt=Aopt3,Zt=C,
            HHt=sigma1opt3,GGt=sigma2,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a6.3$logLik))
}

opt6.3 <- optim(40,func6.3,method = "L-BFGS-B")
initvaropt3<-opt6.3$par*matrix(c(1,0,0,1),nrow=2)

func6.4 <- function(par){
  sigma2opt3<-par*matrix(c(1,0,0,1),nrow=2)
  a6.4<-fkf(a0=initstate,P0=initvaropt3,dt=c(0,0),ct=c(0,0),Tt=Aopt3,Zt=C,
            HHt=sigma1opt3,GGt=sigma2opt3,yt=matrix(c(NOx.diff,NO2.diff),nrow=2))
  return(-(a6.4$logLik))
}

opt6.4 <- optim(25,func6.4,method = "L-BFGS-B",lower=0)
sigma2opt3<-opt6.4$par*matrix(c(1,0,0,1),nrow=2)



aopt3<-fkf(a0=initstate,P0=initvaropt3,dt=c(0,0),ct=c(0,0),Tt=Aopt3,Zt=C,
           HHt=sigma1opt3,GGt=sigma2opt3,
           yt=matrix(c(NOx.diff,NO2.diff),nrow=2))

aopt3$logLik

fcdiffstatesopt3 <- matrix(NA,2,25)
fcsigma1opt3 <- rep(list(matrix(NA,2,2)),25)
for (i in 1:24){
  fcdiffstatesopt3[,1] <- t(t(aopt3$at[,560]))
  fcdiffstatesopt3[,i+1] <- Aopt3%*%fcdiffstatesopt3[,i]
}
for (i in 1:24){
  fcsigma1opt3[[1]] <- aopt3$Pt[,,560]
  fcsigma1opt3[[i+1]] <- Aopt3%*%fcsigma1opt3[[i]]%*%t(Aopt3)+sigma1opt3
}      

fcdiffstatesopt3<-fcdiffstatesopt3[,2:25]
fcsigma1opt3 <- fcsigma1opt3[2:25]
fcdiffopt3 <-C%*%fcdiffstatesopt3
fcsigma2opt3 <- rep(list(matrix(NA,2,2)),25)
for (i in 1:24){
  fcsigma2opt3[[1]] <- fcsigma1opt3[[1]]
  fcsigma2opt3[[i+1]] <- C%*%fcsigma1opt3[[i]]%*%t(C)+sigma2opt3
} 
fcsigma2opt3 <- fcsigma2opt3[2:25]

fcNOxopt3 <- numeric(24)
fcNO2opt3 <- numeric(24)
predintdiffNOxopt3.lower <- numeric(24)
predintdiffNOxopt3.upper <- numeric(24) 
predintdiffNO2opt3.lower <- numeric(24)
predintdiffNO2opt3.upper <- numeric(24)

for (i in 1:24){
  predintdiffNOxopt3.lower[i] <- fcdiffopt3[1,i]-qt(0.975,597)*
    sqrt(aopt3$Pt[1,1,560]+fcsigma2opt3[[i]][1,1])
  predintdiffNOxopt3.upper[i] <- fcdiffopt3[1,i]+qt(0.975,597)*
    sqrt(aopt3$Pt[1,1,560]+fcsigma2opt3[[i]][1,1])
  predintdiffNO2opt3.lower[i] <- fcdiffopt3[2,i]-qt(0.975,597)*
    sqrt(aopt3$Pt[2,2,560]+fcsigma2opt3[[i]][2,2])
  predintdiffNO2opt3.upper[i] <- fcdiffopt3[2,i]+qt(0.975,597)*
    sqrt(aopt3$Pt[2,2,560]+fcsigma2opt3[[i]][2,2])
}

predintNOxopt3.lower <- numeric(24)
predintNOxopt3.upper <- numeric(24)
predintNO2opt3.lower <- numeric(24)
predintNO2opt3.upper <- numeric(24)

for (i in 560:583){
  fcNOxopt3[i-559] <- fcdiffopt3[1,i-559]+train.data[i,3]
  fcNO2opt3[i-559] <- fcdiffopt3[2,i-559]+train.data[i,2]
  predintNOxopt3.lower[i-559] <- predintdiffNOxopt3.lower[i-559]+train.data[i,3]
  predintNOxopt3.upper[i-559] <- predintdiffNOxopt3.upper[i-559]+train.data[i,3]
  predintNO2opt3.lower[i-559] <- predintdiffNO2opt3.lower[i-559]+train.data[i,2]
  predintNO2opt3.upper[i-559] <- predintdiffNO2opt3.upper[i-559]+train.data[i,2]
}

plot(test.data[,3],type="l",col="red",xlim=c(1,24),ylim=c(-95.8874,423.7312),
     main="Predictions for the final 24 hours for NOx",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(fcNOxopt3,col="blue")
lines(predintNOxopt3.lower,col="green",lty=2)
lines(predintNOxopt3.upper,col="green",lty=2)
legend("topleft",legend=c("95% prediction interval for optimised predictions",
                          "Optimised predictions","Testing data"),
       col=c("green","blue","red"),lty=c(2,1,1),cex=0.8)

plot(test.data[,2],type="l",col="red",xlim=c(1,24),ylim=c(-89.55237,250),
     main="Predictions for the final 24 hours for NO2",
     xlab="Observation",ylab="Concentration (micrograms/cubic meter)")
lines(fcNO2opt3,col="blue")
lines(predintNO2opt3.lower,col="green",lty=2)
lines(predintNO2opt3.upper,col="green",lty=2)
legend("topleft",legend=c("95% prediction interval for optimised predictions",
                          "Optimised predictions","Testing data"),
       col=c("green","blue","red"),lty=c(2,1,1),cex=0.8)

predict.table.opt3.a <- data.frame(test.data[,3],fcNOxopt3,predintNOxopt3.lower,
                                   predintNOxopt3.upper)
names(predict.table.opt3.a) <- c("Actual NOx concentrations",
                                 "NOx predicted concentrations",
                                 "NOx prediction interval lower bound",
                                 "NOx prediction interval upper bound")
predict.table.opt3.a <- predict.table.opt3.a[c(1,2,6,24),]

predict.table.opt3.b<-data.frame(test.data[,2],fcNO2opt3,predintNO2opt3.lower,
                                 predintNO2opt3.upper)
names(predict.table.opt3.b)<-c("Actual NO2 concentrations",
                               "NO2 predicted concentrations",
                               "NO2 prediction interval lower bound",
                               "NO2 prediction interval upper bound")
predict.table.opt3.b <- predict.table.opt3.b[c(1,2,6,24),]

