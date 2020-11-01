## Assignment 3

## Question 3.1
data <- data.frame(read.csv("A3_hcab.csv"))
data[,5] <- c(1:751)
train.data<-data[1:703,]
test.data<-data[704:751,]
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(train.data[,2],type="l",col="red",xlim=c(1,751),ylim=c(1.26,404.53),
     main="Concentrations of the three gasses",xlab="Observation",
     ylab="Concentrations (micrograms/cubic meter)")
lines(train.data[,3],col="green")
lines(train.data[,4],col="blue")
lines(data[704:751,2]~data[704:751,5],col="red",lty=2)
lines(data[704:751,3]~data[704:751,5],col="green",lty=2)
lines(data[704:751,4]~data[704:751,5],col="blue",lty=2)
legend("topleft",legend=c("NO2","NOx","O3"),col=c("red","green","blue"),
       lty=1:1, cex=0.8)
legend("top",legend=c("Training data","Testing data"),col=c("black",
                                                             "black"),
                      lty=c(1,2),cex=0.8)

## Question 3.2

par(mfrow=c(2,1), mgp=c(2,0.7,0), mar=c(3,3,1,1))
acf(train.data[,2],lag.max=96)
pacf(train.data[,2],lag.max=96)

acf(train.data[,3],lag.max=96)
pacf(train.data[,3],lag.max=96)

acf(train.data[,4],lag.max=96)
pacf(train.data[,4],lag.max=96)

## All stationary as ACF decreases exponentially to 0
## All seasonal over 24 hour periods

## transformed ACF and PACF
acf(log(train.data[,2]),lag.max=96)
pacf(log(train.data[,2]),lag.max=96)

acf(log(train.data[,3]),lag.max=96)
pacf(log(train.data[,3]),lag.max=96)

acf(sqrt(train.data[,4]),lag.max=96)
pacf(sqrt(train.data[,4]),lag.max=96)

## Question 3.4
library(car)
library(Hmisc)
par(mfrow=c(1,1),mgp=c(2, 0.7,0), mar=c(3,3,2,1))
m1 <- arima(train.data[,2],order=c(1,0,0),
            seasonal=list(order=c(0,1,1),period=24))
## start at (3,0,0)x(0,1,0) but realise 2nd and 3rd AR variables are
## not neccesary. SMA variable needed
pacf(residuals(m1))
tsdiag(m1)
m1
hist(residuals(m1))
qqPlot(residuals(m1))
cpgram(residuals(m1))
## heavy tails so transformation

m1s <- arima(sqrt(train.data[,2]),order= c(1,0,0),
             seasonal=list(order=c(0,1,1),period=24))
pacf(residuals(m1s))
tsdiag(m1s)
m1s
hist(residuals(m1s))
qqPlot(residuals(m1s))
cpgram(residuals(m1s))

m1l <- arima(log(train.data[,2]),order= c(1,0,0),
             seasonal=list(order=c(0,1,1),period=24))
pacf(residuals(m1l))
tsdiag(m1l)
m1l
hist(residuals(m1l))
qqPlot(residuals(m1l))
cpgram(residuals(m1l))

m1.residuals <- rep(NA,703)
for (i in 1:703){
  if (m1$residuals[i]<0){
    m1.residuals[i]=-1}
  else 
    m1.residuals[i]=1
}
sumofsigns <- 0
for (i in 1:702){
  if (m1.residuals[i]!=m1.residuals[i+1]){
    sumofsigns <- sumofsigns+1}
  else
    sumofsigns <- sumofsigns
}
sumofsigns/703
m1s.residuals <- rep(NA,703)
for (i in 1:703){
  if (m1s$residuals[i]<0){
    m1s.residuals[i]=-1}
  else 
    m1s.residuals[i]=1
}
sumofsigns <- 0
for (i in 1:702){
  if (m1s.residuals[i]!=m1s.residuals[i+1]){
    sumofsigns <- sumofsigns+1}
  else
    sumofsigns <- sumofsigns
}
sumofsigns/703
m1l.residuals <- rep(NA,703)
for (i in 1:703){
  if (m1l$residuals[i]<0){
  m1l.residuals[i]=-1}
  else 
  m1l.residuals[i]=1
}
sumofsigns <- 0
for (i in 1:702){
  if (m1l.residuals[i]!=m1l.residuals[i+1]){
    sumofsigns <- sumofsigns+1}
  else
    sumofsigns <- sumofsigns
}
sumofsigns/703
## log has a lower AIC and better ratio of sign changes, closer to 1/2
## squareroot has better qqplot
binconf(331,703)
## within 95% confidence interval so can be considered white noise


m2 <- arima(train.data[,3],order= c(1,0,0),
            seasonal=list(order=c(1,1,0),period=24))
## start at (3,0,0)x(0,1,0) but realise 2nd and 3rd AR variables are
## not neccesary. SAR variable needed
pacf(residuals(m2))
tsdiag(m2)
m2
hist(residuals(m2))
qqPlot(residuals(m2))
cpgram(residuals(m2))
## heavy tails so transformation needed

m2s <- arima(sqrt(train.data[,3]),order= c(1,0,0),
             seasonal=list(order=c(1,1,0),period=24))
pacf(residuals(m2s))
tsdiag(m2s)
m2s
hist(residuals(m2s))
qqPlot(residuals(m2s))
cpgram(residuals(m2s))

m2l <- arima(log(train.data[,3]),order= c(2,0,0),
             seasonal=list(order=c(0,1,1),period=24))
pacf(residuals(m2l))
tsdiag(m2l)
m2l
hist(residuals(m2l))
qqPlot(residuals(m2l))
cpgram(residuals(m2l))

m2.residuals <- rep(NA,703)
for (i in 1:703){
  if (m2$residuals[i]<0){
    m2.residuals[i]=-1}
  else 
    m2.residuals[i]=1
}
sumofsigns <- 0
for (i in 1:702){
  if (m2.residuals[i]!=m2.residuals[i+1]){
    sumofsigns <- sumofsigns+1}
  else
    sumofsigns <- sumofsigns
}
sumofsigns/703
m2s.residuals <- rep(NA,703)
for (i in 1:703){
  if (m2s$residuals[i]<0){
    m2s.residuals[i]=-1}
  else 
    m2s.residuals[i]=1
}
sumofsigns <- 0
for (i in 1:702){
  if (m2s.residuals[i]!=m2s.residuals[i+1]){
    sumofsigns <- sumofsigns+1}
  else
    sumofsigns <- sumofsigns
}
sumofsigns/703
m2l.residuals <- rep(NA,703)
for (i in 1:703){
  if (m2l$residuals[i]<0){
    m2l.residuals[i]=-1}
  else 
    m2l.residuals[i]=1
}
sumofsigns <- 0
for (i in 1:702){
  if (m2l.residuals[i]!=m2l.residuals[i+1]){
    sumofsigns <- sumofsigns+1}
  else
    sumofsigns <- sumofsigns
}
sumofsigns/703
## both have same ratio of sign changes, close to 1/2
## log has a lower AIC and better qqplot
binconf(349,703)
## within 95% confidence interval so can be considered white noise


m3 <- arima(train.data[,4],order= c(2,0,0),
            seasonal=list(order=c(0,1,1),period=24))
## start at (3,0,0)x(0,1,0) but realise 3rd AR variable not neccesary.
## SMA variable needed
pacf(residuals(m3))
tsdiag(m3)
m3
hist(residuals(m3))
qqPlot(residuals(m3))
cpgram(residuals(m3))

m3s <- arima(sqrt(train.data[,4]),order= c(2,0,0),
             seasonal=list(order=c(0,1,1),period=24))
pacf(residuals(m3s))
tsdiag(m3s)
m3s
hist(residuals(m3s))
qqPlot(residuals(m3s))
cpgram(residuals(m3s))

m3l <- arima(log(train.data[,4]),order= c(3,0,0),
             seasonal=list(order=c(0,1,1),period=24))
pacf(residuals(m3l))
tsdiag(m3l)
m3l
hist(residuals(m3l))
qqPlot(residuals(m3l))
cpgram(residuals(m3l))

m3.residuals <- rep(NA,703)
for (i in 1:703){
  if (m3$residuals[i]<0){
    m3.residuals[i]=-1}
  else 
    m3.residuals[i]=1
}
sumofsigns <- 0
for (i in 1:702){
  if (m3.residuals[i]!=m3.residuals[i+1]){
    sumofsigns <- sumofsigns+1}
  else
    sumofsigns <- sumofsigns
}
sumofsigns/703
m3s.residuals <- rep(NA,703)
for (i in 1:703){
  if (m3s$residuals[i]<0){
    m3s.residuals[i]=-1}
  else 
    m3s.residuals[i]=1
}
sumofsigns <- 0
for (i in 1:702){
  if (m3s.residuals[i]!=m3s.residuals[i+1]){
    sumofsigns <- sumofsigns+1}
  else
    sumofsigns <- sumofsigns
}
sumofsigns/703
m3l.residuals <- rep(NA,703)
for (i in 1:703){
  if (m3l$residuals[i]<0){
    m3l.residuals[i]=-1}
  else 
    m3l.residuals[i]=1
}
sumofsigns <- 0
for (i in 1:702){
  if (m3l.residuals[i]!=m3l.residuals[i+1]){
    sumofsigns <- sumofsigns+1}
  else
    sumofsigns <- sumofsigns
}
sumofsigns/703
## squareroot has better ratio of sign changes, is closer to 1/2. Also
## has better qqplot
## log has lower AIC
binconf(346,703)
## within 95% confidence interval so can be considered white noise



## Question 3.5
pred1l <- predict(m1l,n.ahead=48)
pred1l.int <- rep(NA,48)
for(i in 1:48){
  pred1l.int[i] <-
    qt(0.975, 702)*(pred1l$se[i])
}
pred1l.upp <- pred1l$pred+pred1l.int
pred1l.low <- pred1l$pred-pred1l.int
plot(test.data[,2],type="l",col="red",xlim=c(1,48),ylim=c(4.6,118.8689),
     main="Predictions for concentration of NO2",xlab="Observation",
     ylab="Concentration (micrograms/cubic meter)")
lines(exp(pred1l$pred)~c(1:48),col="blue",lty=1)
lines(exp(pred1l.upp)~c(1:48),col="green",lty=1)
lines(exp(pred1l.low)~c(1:48),col="green",lty=1)
legend("topleft",legend=c("Testing data","Predictions",
                          "95% prediction interval"),
       col=c("red","blue","green"),lty=1:1, cex=0.8)
tab1l <- data.frame(test.data[,2],exp(pred1l$pred),exp(pred1l.low),
           exp(pred1l.upp))
names(tab1l)[1]<-paste("Testing data")
names(tab1l)[2]<-paste("Model predictions")
names(tab1l)[3]<-paste("95% prediction interval lower bound")
names(tab1l)[4]<-paste("95% prediction interval upper bound")
tab1l <- tab1l[c(1,6,12,24,48),]


pred2l <- predict(m2l,n.ahead=48)
pred2l.int <- rep(NA,48)
for(i in 1:48){
  pred2l.int[i] <-
    qt(0.975, 702)*(pred2l$se[i])
}
pred2l.upp <- pred2l$pred+pred2l.int
pred2l.low <- pred2l$pred-pred2l.int
plot(test.data[,3],type="l",col="red",xlim=c(1,48),ylim=c(5.8,404.53),
     main="Predictions for concentration of NOx",xlab="Observation",
     ylab="Concentration (micrograms/cubic meter)")
lines(exp(pred2l$pred)~c(1:48),col="blue",lty=1)
lines(exp(pred2l.upp)~c(1:48),col="green",lty=1)
lines(exp(pred2l.low)~c(1:48),col="green",lty=1)
legend("topleft",legend=c("Testing data","Predictions",
                          "95% prediction interval"),
       col=c("red","blue","green"),lty=1:1, cex=0.8)
tab2l <- data.frame(test.data[,3],exp(pred2l$pred),exp(pred2l.low),
                    exp(pred2l.upp))
names(tab2l)[1]<-paste("Testing data")
names(tab2l)[2]<-paste("Model predictions")
names(tab2l)[3]<-paste("95% prediction interval lower bound")
names(tab2l)[4]<-paste("95% prediction interval upper bound")
tab2l <- tab2l[c(1,6,12,24,48),]




pred3s <- predict(m3s,n.ahead=48)
pred3s.int <- rep(NA,48)
for(i in 1:48){
  pred3s.int[i] <-
    qt(0.975, 702)*(pred3s$se[i])
}
pred3s.upp <- pred3s$pred+pred3s.int
pred3s.low <- pred3s$pred-pred3s.int
plot(test.data[,4],type="l",col="red",xlim=c(1,48),ylim=c(-10,89.12664),
     main="Predictions for concentration of O3",xlab="Observation",
     ylab="Concentration (micrograms/cubic meter)")
lines(((pred3s$pred)^2)~c(1:48),col="blue",lty=1)
lines(((pred3s.upp)^2)~c(1:48),col="green",lty=1)
lines(((pred3s.low)^2)~c(1:48),col="green",lty=1)
legend("bottomleft",legend=c("Testing data","Predictions",
                          "95% prediction interval"),
       col=c("red","blue","green"),lty=1:1, cex=0.8)
tab3s <- data.frame(test.data[,4],(pred3s$pred)^2,(pred3s.low)^2,
                    (pred3s.upp)^2)
names(tab3s)[1]<-paste("Testing data")
names(tab3s)[2]<-paste("Model predictions")
names(tab3s)[3]<-paste("95% prediction interval lower bound")
names(tab3s)[4]<-paste("95% prediction interval upper bound")
tab3s <- tab3s[c(1,6,12,24,48),]


## Question 3.6
multdata <- cbind(train.data[,2],train.data[,3],train.data[,4])
acf(multdata,lag.max=96)
multdatatr <- cbind(log(train.data[,2]),log(train.data[,3]),
                    sqrt(train.data[,4]))
acf(multdatatr,lag.max=96)


## Question 3.7
library(marima)

dd <- define.dif(multdata, difference = c(1,24, 2,24, 3,24))
Mod <- define.model(kvar=3,ar=c(1),ma=c(0))
Marima <- marima(dd$y.dif,ar.pattern=Mod$ar.pattern,
                  ma.pattern=Mod$ma.pattern, Plot='log.det')
short.form(Marima$ar.estimates, leading=FALSE)
short.form(Marima$ma.estimates, leading=FALSE)
Marima$ar.fvalues[,,2]
Marima$ma.fvalues[,,1]

source("step.slow.p.marima_2017.R")
ss <- step.slow.p(object = Marima, data = dd$y.dif)

short.form(ss$ar.estimates)
short.form(ss$ar.fvalues)
short.form(ss$ar.stdv) 


multdatatr <- cbind(log(train.data[,2]),log(train.data[,3]),
                    sqrt(train.data[,4]))

ddtr <- define.dif(multdatatr, difference = c(1,24, 2,24, 3,24))
Mod2 <- define.model(kvar=3,ar=c(1),ma=c(0))
Marima2 <- marima(ddtr$y.dif,ar.pattern=Mod2$ar.pattern,
                  ma.pattern=Mod2$ma.pattern, Plot='log.det')
short.form(Marima2$ar.estimates, leading=FALSE)
short.form(Marima2$ma.estimates, leading=FALSE)
Marima2$ar.fvalues[,,2]
Marima2$ma.fvalues[,,1]

ss2 <- step.slow.p(object = Marima2, data = ddtr$y.dif)

short.form(ss2$ar.estimates)
short.form(ss2$ar.fvalues)
short.form(ss2$ar.stdv) 

## Question 3.8

predtr <- arma.forecast(series = rbind(multdatatr[1:703,], matrix(NA, nrow=48, ncol=3)), 
                      marima = ss2, nstart = 703, nstep = 48, dif.poly = dd$dif.poly)

fcast2 <- t(predtr$forecasts)[704:751,]

plot(test.data[,2], type="l",col="red", ylim=c(2.423068,136.5793))
lines(exp(fcast2[,1]),col="blue")

sd1tr <- sqrt(predtr$pred.var[1,1,])
lims1tr.lower <- (fcast2[,1]-sd1tr*qt(p=0.975,702))
lims1tr.upper <- (fcast2[,1]+sd1tr*qt(p=0.975,702))
lines(exp(lims1tr.lower),col="green")
lines(exp(lims1tr.upper),col="green")
legend("topleft",legend=c("Testing data","Predictions",
                          "95% prediction interval"),
       col=c("red","blue","green"),lty=1:1, cex=0.8)
tab1mltr <- data.frame(test.data[,2],exp(fcast2[,1]),exp(lims1tr.lower),
                    exp(lims1tr.upper))
names(tab1mltr)[1]<-paste("Testing data")
names(tab1mltr)[2]<-paste("Model predictions")
names(tab1mltr)[3]<-paste("95% prediction interval lower bound")
names(tab1mltr)[4]<-paste("95% prediction interval upper bound")
tab1mltr <- tab1mltr[c(1,6,12,24,48),]


plot(test.data[,3], type="l",col="red", ylim=c(3.821881,404.53))
lines(exp(fcast2[,2]),col="blue")

sd2tr <- sqrt(predtr$pred.var[2,2,])
lims2tr.lower <- (fcast2[,2]-sd2tr*qt(p=0.975,702))
lims2tr.upper <- (fcast2[,2]+sd2tr*qt(p=0.975,702))
lines(exp(lims2tr.lower),col="green")
lines(exp(lims2tr.upper),col="green")
legend("topleft",legend=c("Testing data","Predictions",
                          "95% prediction interval"),
       col=c("red","blue","green"),lty=1:1, cex=0.8)
tab2mltr <- data.frame(test.data[,3],exp(fcast2[,2]),exp(lims2tr.lower),
                       exp(lims2tr.upper))
names(tab2mltr)[1]<-paste("Testing data")
names(tab2mltr)[2]<-paste("Model predictions")
names(tab2mltr)[3]<-paste("95% prediction interval lower bound")
names(tab2mltr)[4]<-paste("95% prediction interval upper bound")
tab2mltr <- tab2mltr[c(1,6,12,24,48),]


plot(test.data[,4], type="l",col="red", ylim=c(0.6414558,140.6499))
lines((fcast2[,3])^2,col="blue")

sd3tr <- sqrt(predtr$pred.var[3,3,])
lims3tr.lower <- (fcast2[,3]-sd3tr*qt(p=0.975,702))
lims3tr.upper <- (fcast2[,3]+sd3tr*qt(p=0.975,702))
lines((lims3tr.lower)^2,col="green")
lines((lims3tr.upper)^2,col="green")
legend("topleft",legend=c("Testing data","Predictions",
                          "95% prediction interval"),
       col=c("red","blue","green"),lty=1:1, cex=0.8)
tab3mltr <- data.frame(test.data[,4],(fcast2[,3]^2),(lims3tr.lower^2),
                       (lims3tr.upper^2))
names(tab3mltr)[1]<-paste("Testing data")
names(tab3mltr)[2]<-paste("Model predictions")
names(tab3mltr)[3]<-paste("95% prediction interval lower bound")
names(tab3mltr)[4]<-paste("95% prediction interval upper bound")
tab3mltr <- tab3mltr[c(1,6,12,24,48),]
