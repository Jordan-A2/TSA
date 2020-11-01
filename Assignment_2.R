## Assignment 2

## Question 2.1
## 3)
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
sim <- matrix(NA,200,10)
for (i in 1:10){
  sim[,i] <- arima.sim(model=list(ar=0.9,ma=c(2,0.5),order=c(1,0,2)),
                      n =200)
  plot(sim[,i],type="l",col="red",
       main=paste("Arima(1,0,2) simulation",i),xlab="Time",ylab="X_t")
}


matplot(sim,type="o",lty=1,pch=1,col=rainbow(11),
        main="ARIMA(1,0,2) simulations",xlab="Time",
     ylab="X_t",cex=0.4,ylim=c(min(sim),max(sim)))


## 4)
lmax <- 23 ## lag.max
acf.dat <- matrix(NA,ncol=10,nrow=(lmax+1))
for (i in 1:10){
  acf.dat[,i]<-acf(sim[,i],lag.max=lmax)$acf
}

matplot(acf.dat,type="o",lty=1,pch=1,col=rainbow(11),xlab="Lags",
        ylab="Auto-correlation",main="Auto-correlations for the 10 simulations")

## 5)
lmax <- 23 ## lag.max
pacf.dat <- matrix(NA,ncol=10,nrow=(lmax+1))
for (i in 1:10){
  pacf.dat[,i]<-pacf(sim[,i],lag.max=lmax+1)$acf
}

matplot(pacf.dat,type="o",lty=1,pch=1,col=rainbow(11),xlab="Lags",
        ylab="Partial auto-correlation",
        main="Partial auto-correlations for the 10 simulations")

## 6)
for (i in 1:10){
  print(var(sim[,i]))
}


## Question 2.2
Time <- c(1:16)
Cons <- c(190,208,213,223,237,214,221,201,191,184,184,189,188,207,221,225)
Pred.time <- c(17,18)
Pred.cons <- c(232.6,217.3)
Pred.int.lower <- c(232.6-(qnorm(0.975)*7),
                    217.3-(qnorm(0.975)*7*sqrt(1.25)))
Pred.int.upper <- c(232.6+(qnorm(0.975)*7),
                     217.3+(qnorm(0.975)*7*sqrt(1.25)))
plot(Cons~Time,type="o",pch=19,cex=0.5,main="Electricity consumption for Copenhagen",
     xlab="Date",ylab="Electricity consumption (GWh)",col="red",
     xlim=c(1,18),ylim=c(184,246))
points(Pred.cons~Pred.time,pch=19,cex=0.5,col="blue")
lines(Pred.cons~Pred.time,col="blue")
lines(Pred.int.lower~Pred.time,col="green")
lines(Pred.int.upper~Pred.time,col="green")
legend("top",legend=c("Consumption data","Consumption predictions",
                          "95% prediction interval"),
       col=c("red","blue","green"),lty=1:1, cex=0.8)

predict.table <- data.frame(Pred.time,Pred.cons,Pred.int.lower,
                            Pred.int.upper)



## Question 2.3
## 1)
Mod_1 <- arima.sim(model=list(ar=-0.85),n = 1000)
plot(Mod_1,col="red",main="ARIMA(1,0,0)X(0,0,0) simulation",
     xlab="Time",ylab="X_t",type="l")

acf(Mod_1)

pacf(Mod_1)

## 2)
Mod_2 <- arima.sim(model=list(ar=c(0,0,0,0,0,0,0,0,0,0,0,0.8)),n=1000)
plot(Mod_2,col="red",main="ARIMA(0,0,0)X(1,0,0) simulation",
     xlab="Time",ylab="X_t",type="l")

acf(Mod_2)

pacf(Mod_2)


## 3)
Mod_3 <- arima.sim(model=list(ar=0.8,
                              ma=c(0,0,0,0,0,0,0,0,0,0,0,0.7)),n=1000)
plot(Mod_3,col="red",main="ARIMA(1,0,0)X(0,0,1) simulation",
     xlab="Time",ylab="X_t",type="l")

acf(Mod_3)

pacf(Mod_3)


## 4)
Mod_4 <- arima.sim(model=list(ar=c(0.7,0,0,0,0,0,0,0,0,0,0,-0.8,0.7*0.8)),
                              n =1000)
plot(Mod_4,col="red",main="ARIMA(1,0,0)X(1,0,0) simulation",
     xlab="Time",ylab="X_t",type="l")

acf(Mod_4)

pacf(Mod_4)


## 5)
Mod_5 <- arima.sim(model=list(ar=c(-0.6,-0.6,0,0,0,0,0,0,0,0,0,-0.8,-(0.6*0.8),-(0.6*0.8))),
                   n=1000)
plot(Mod_5,col="red",main="ARIMA(2,0,0)X(1,0,0) simulation",
     xlab="Time",ylab="X_t",type="l")

acf(Mod_5)

pacf(Mod_5)


## 6)
Mod_6 <- arima.sim(model=list(ma=c(-0.7,0,0,0,0,0,0,0,0,0,0,0.8,-(0.7*0.8))),
                   n=1000)
plot(Mod_6,col="red",main="ARIMA(0,0,1)X(0,0,1) simulation",
     xlab="Time",ylab="X_t",type="l")

acf(Mod_6)

pacf(Mod_6)

