## Question 1
df <- read.table("A1_annual.txt", header = TRUE)
df$sh <- NULL

## defining the training and testing data
Training_data <- df[-c(165,166,167,168,169), ]
Testing_data <- df[c(165,166,167,168,169), ]

par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(Training_data$nh~Training_data$year,
     main="Temperature anomalies as a function of time",xlab="Year",
     ylab="Temperature anomality",pch=19,col="red",cex=0.5,
     xlim=c(1850,2018),ylim=c(-0.668,1.064))
points(Testing_data$nh~Testing_data$year,pch=19,col="blue",cex=0.5)
legend("topleft",legend=c("Training data","Testing data"),
       col=c("red","blue"),lty=1:1, cex=0.8)


## Question 2

## initialising for the first 3 values of temperature anomality
first <- 3
pred.nh <- numeric(164)
pred.nh[1:3] <- NA
first_data <- c(Training_data$nh[1:3])
time <- Training_data$year-1849
L <- matrix(c(1,1,0,1),ncol=2)
LInv <- solve(L)
X3 <- cbind(1,(time[1:3])-first)
F_3 <- t(X3)%*%X3
h_3 <- t(X3)%*%first_data
theta3 <- solve(F_3)%*%h_3

##  creating space
theta <- matrix(0,2,164)
F <- rep(list(matrix(0,2,2)),164)
h <- matrix(0,2,164)

## one-step predicitions for the global linear trend model
for (j in (first):163){
  theta[,first] <- theta3
  F[[first]] <- F_3
  h[,first] <- h_3
  pred.nh[j+1] <-  c(1,1)%*%theta[,j]
  F[[j+1]] <- F[[j]]+(c(1,-j))%*%t(c(1,-j))
  h[,j+1] <- LInv%*%h[,j]+(c(1,0))*Training_data$nh[j+1]
  theta[,j+1] <- solve(F[[j+1]])%*%h[,j+1]
}

plot(pred.nh~Training_data$year, type="l",main="One-step predictions",
     xlab="Year",ylab="Temperature anomality",
     col="red",xlim=c(1850,2018),ylim=c(-0.668,1.064))
lines(Training_data$nh~Training_data$year,col="blue")
legend("topleft",legend=c("Training data","One-step predictions"),
       col=c("blue","red"),lty=1:1, cex=0.8)

## calculating the prediction errors
pred.err <- Training_data$nh-pred.nh  
pred.err[1:3] <- NA
plot(pred.err~Training_data$year,type="l",
     main="One-step prediction errors",xlab="Year",
     ylab="Prediction error", col="red",xlim=c(1850,2018),
     ylim=c(-0.668,1.064))

## calculating the predictions of the testing data
pred.test <- c(c(1,1)%*%theta[,164],c(1,2)%*%theta[,164],
               c(1,3)%*%theta[,164],c(1,4)%*%theta[,164],
               c(1,5)%*%theta[,164])

## calculating the variance of the prediction error for the 164th value
## where values 3:164 have been used
X <- cbind(1,-161:0)
S <- (Training_data$nh[3:164]-X%*%theta[,164])
var <- (t(S)%*%S) / (160)


## 95% prediction intervals for the testing data
predict.int <- numeric(5)
for(i in 1:5){
  predict.int[i] <-
    qt(0.975, 160)*sqrt(var)* 
    sqrt((1+t(c(1,i))%*%solve(F[[164]])%*%c(1,i)))
}

pred.test.lower <- pred.test - predict.int
pred.test.upper <- pred.test + predict.int

plot(pred.nh~Training_data$year, type="l",
     main="One-step predictions with testing data predictions",
     xlab="Year",ylab="Temperature anomality",
     col="red",xlim=c(1850,2018),ylim=c(-0.668,1.064))
lines(Training_data$nh~Training_data$year,col="black")
lines(pred.test~Testing_data$year,col="blue")
lines(pred.test.lower~Testing_data$year,col="green")
lines(pred.test.upper~Testing_data$year,col="green")
points(Testing_data$nh~Testing_data$year,pch=19,cex=0.5,col="orange")
legend("topleft",legend=c("Training data","One-step predictions",
                          "Testing data predictions",
                          "95% prediction interval","Testing data"),
       col=c("black","red","blue","green","orange"),lty=1:1, cex=0.8)


## table comparing the global linear trend model predictions of the
## testing data and the actual testing data. Also includes the upper
## and lower bounds from the prediction interval
predict.table <- data.frame(Testing_data$year,pred.test,
                            pred.test.lower,pred.test.upper
                            ,Testing_data$nh)
names(predict.table) <- c("Year","Testing data predictions",
                          "95% prediction interval lower bound",
                          "95% prediction interval upper bound",
                          "Testing data")


## Question 3
## initialising for the first 3 values of temperature anomality
first <- 3
Lpred.nh <- numeric(164)
Lpred.nh[1:3] <- NA
first_data <- c(Training_data$nh[1:3])
time <- Training_data$year-1849
L <- matrix(c(1,1,0,1),ncol=2)
LInv <- solve(L)

## creating the initial epsilon matrix
lambda_init <- matrix(c((1/0.8^2),0,0, 0,(1/0.8),0, 0,0,1),nrow=3)

X3 <- cbind(1,(time[1:3])-first)
LF_3 <- t(X3)%*%solve(lambda_init)%*%X3
Lh_3 <- t(X3)%*%solve(lambda_init)%*%first_data
Ltheta3 <- solve(LF_3)%*%Lh_3

##creating space
Ltheta <- matrix(0,2,164)
LF <- rep(list(matrix(0,2,2)),164)
Lh <- matrix(0,2,164)

## one-step predicitions for the local linear trend model
for (j in (first):163){
  Ltheta[,first] <- Ltheta3
  LF[[first]] <- LF_3
  Lh[,first] <- Lh_3
  Lpred.nh[j+1] <-  c(1,1)%*%Ltheta[,j]
  LF[[j+1]] <- LF[[j]]+(0.8^j)*(c(1,-j))%*%t(c(1,-j))
  Lh[,j+1] <- (0.8)*LInv%*%Lh[,j]+(c(1,0))*Training_data$nh[j+1]
  Ltheta[,j+1] <- solve(LF[[j+1]])%*%Lh[,j+1]
}

plot(Lpred.nh~Training_data$year, type="l",
     main="Local one-step predictions",
     xlab="Year",ylab="Temperature anomality",
     col="red",xlim=c(1850,2018),ylim=c(-0.668,1.064))
lines(Training_data$nh~Training_data$year,col="blue")
legend("topleft",legend=c("Training data","Local one-step predictions"),
       col=c("blue","red"),lty=1:1, cex=0.8)

## calculating the prediction errors
Lpred.err <- Training_data$nh-Lpred.nh  
Lpred.err[1:3] <- NA
plot(Lpred.err~Training_data$year,type="l",
     main="Local one-step prediction errors",xlab="Year",
     ylab="Prediction error", col="red",xlim=c(1850,2018),
     ylim=c(-0.668,1.064))

## calculating the predictions of the testing data
Lpred.test <- c(c(1,1)%*%Ltheta[,164],c(1,2)%*%Ltheta[,164],
               c(1,3)%*%Ltheta[,164],c(1,4)%*%Ltheta[,164],
               c(1,5)%*%Ltheta[,164])

## calculating the variance of the prediction error for the 164th value
## where values 3:164 have been used. The divisor in Lvar is calculated
## using T-p where T=(1/1-0.8) and p=2
X <- cbind(1,-161:0)
LS <- (Training_data$nh[3:164]-X%*%Ltheta[,164])
Lepsilon <- matrix(0,162,162)
for (i in 162:1){
  Lepsilon[163-i,163-i] <- 1/(0.8^(i-1))
}
Lvar <- (t(LS)%*%solve(Lepsilon)%*%LS) / (3)

## 95% prediction intervals for the testing data
Lpredict.int <- numeric(5)
for(i in 1:5){
  Lpredict.int[i] <-
    qt(0.975, 160)*sqrt(Lvar)* 
    sqrt((1+t(c(1,i))%*%solve(LF[[164]])%*%c(1,i)))
}

Lpred.test.lower <- Lpred.test - Lpredict.int
Lpred.test.upper <- Lpred.test + Lpredict.int

plot(Lpred.nh~Training_data$year, type="l",
     main="Local one-step predictions with testing data predictions",
     xlab="Year",ylab="Temperature anomality",
     col="red",xlim=c(1850,2018),ylim=c(-0.668,1.064))
lines(Training_data$nh~Training_data$year,col="black")
lines(Lpred.test~Testing_data$year,col="blue")
lines(Lpred.test.lower~Testing_data$year,col="green")
lines(Lpred.test.upper~Testing_data$year,col="green")
points(Testing_data$nh~Testing_data$year,pch=19,cex=0.5,col="orange")
legend("topleft",legend=c("Training data","Local one-step predictions",
                          "Testing data predictions",
                          "95% prediction interval","Testing data"),
       col=c("black","red","blue","green","orange"),lty=1:1, cex=0.8)


## table comparing the local linear trend model predictions of the
## testing data and the actual testing data. Also includes the upper
## and lower bounds from the prediction interval
Lpredict.table <- data.frame(Testing_data$year,Lpred.test,
                            Lpred.test.lower,Lpred.test.upper
                            ,Testing_data$nh)
names(Lpredict.table) <- c("Year","Testing data predictions",
                          "95% prediction interval lower bound",
                          "95% prediction interval upper bound",
                          "Testing data")


## Question 4
## minimising the sum of squared errors with values of lambda between
## 0.5 and 1 by intervals of 0.01
## running the same code as question 1.3 but with an extra for loop
## around the code where lambda changes value

## creating space for the sum of square errors for each lambda
sum_sq_err <- numeric(51)
for (k in seq(0.5,1,0.01)){
first <- 3
OLpred.nh <- numeric(164)
OLpred.nh[1:3] <- NA
first_data <- c(Training_data$nh[1:3])
time <- Training_data$year-1849
L <- matrix(c(1,1,0,1),ncol=2)
LInv <- solve(L)
Olambda_init <- matrix(c((1/k^2),0,0, 0,(1/k),0, 0,0,1),nrow=3)
X3 <- cbind(1,(time[1:3])-first)
OLF_3 <- t(X3)%*%solve(Olambda_init)%*%X3
OLh_3 <- t(X3)%*%solve(Olambda_init)%*%first_data
OLtheta3 <- solve(OLF_3)%*%OLh_3
OLtheta <- matrix(0,2,164)
OLF <- rep(list(matrix(0,2,2)),164)
OLh <- matrix(0,2,164)
for (j in (first):163){
  OLtheta[,first] <- OLtheta3
  OLF[[first]] <- OLF_3
  OLh[,first] <- OLh_3
  OLpred.nh[j+1] <-  c(1,1)%*%OLtheta[,j]
  OLF[[j+1]] <- OLF[[j]]+(k^j)*(c(1,-j))%*%t(c(1,-j))
  OLh[,j+1] <- (k)*LInv%*%OLh[,j]+(c(1,0))*Training_data$nh[j+1]
  OLtheta[,j+1] <- solve(OLF[[j+1]])%*%OLh[,j+1]
}
## not using the first 3 values as they are not one-step predictions
## and not using the first 5 one-step predictions as they are used
## as the burn in period
OLpred.err <- (Training_data$nh[9:164]-OLpred.nh[9:164])^2
sum_sq_err[100*k-49] <- sum(OLpred.err)
}
## plotting the data
## gives a value of 0 at about 0.58 but this is an error (not with the
## code but i think with the processing capabilities of my computer)
k <- c(seq(0.5,1,0.01))
plot(sum_sq_err~k,type="l",col="red",main="Sum of squared errors for values of lambda",
     xlab="Lambda",ylab="sum of squared errors")


## running the same code over a new range to get a more precise value
## for the optimal choice of lambda

## creating space for the new sum of square errors for each lambda
sum_sq_err <- numeric(16)
for (k in seq(0.8,0.95,0.01)){
  first <- 3
  OLpred.nh <- numeric(164)
  OLpred.nh[1:3] <- NA
  first_data <- c(Training_data$nh[1:3])
  time <- Training_data$year-1849
  L <- matrix(c(1,1,0,1),ncol=2)
  LInv <- solve(L)
  Olambda_init <- matrix(c((1/k^2),0,0, 0,(1/k),0, 0,0,1),nrow=3)
  X3 <- cbind(1,(time[1:3])-first)
  OLF_3 <- t(X3)%*%solve(Olambda_init)%*%X3
  OLh_3 <- t(X3)%*%solve(Olambda_init)%*%first_data
  OLtheta3 <- solve(OLF_3)%*%OLh_3
  OLtheta <- matrix(0,2,164)
  OLF <- rep(list(matrix(0,2,2)),164)
  OLh <- matrix(0,2,164)
  for (j in (first):163){
    OLtheta[,first] <- OLtheta3
    OLF[[first]] <- OLF_3
    OLh[,first] <- OLh_3
    OLpred.nh[j+1] <-  c(1,1)%*%OLtheta[,j]
    OLF[[j+1]] <- OLF[[j]]+(k^j)*(c(1,-j))%*%t(c(1,-j))
    OLh[,j+1] <- (k)*LInv%*%OLh[,j]+(c(1,0))*Training_data$nh[j+1]
    OLtheta[,j+1] <- solve(OLF[[j+1]])%*%OLh[,j+1]
  }
  ## not using the first 3 values as they are not one-step predictions
  ## and not using the first 5 one-step predictions as they are used
  ## as the burn in period
  OLpred.err <- (Training_data$nh[9:164]-OLpred.nh[9:164])^2
  sum_sq_err[100*k-79] <- sum(OLpred.err)
}
k <- c(seq(0.8,0.95,0.01))
plot(sum_sq_err~k,type="l",col="red",main="Sum of squared errors for values of lambda",
     xlab="Lambda",ylab="sum of squared errors")


## running the same code over a new range to get a more precise value
## for the optimal choice of lambda

## creating space for the new sum of square errors for each lambda
sum_sq_err <- numeric(11)
for (k in seq(0.845,0.855,0.001)){
  first <- 3
  OLpred.nh <- numeric(164)
  OLpred.nh[1:3] <- NA
  first_data <- c(Training_data$nh[1:3])
  time <- Training_data$year-1849
  L <- matrix(c(1,1,0,1),ncol=2)
  LInv <- solve(L)
  Olambda_init <- matrix(c((1/k^2),0,0, 0,(1/k),0, 0,0,1),nrow=3)
  X3 <- cbind(1,(time[1:3])-first)
  OLF_3 <- t(X3)%*%solve(Olambda_init)%*%X3
  OLh_3 <- t(X3)%*%solve(Olambda_init)%*%first_data
  OLtheta3 <- solve(OLF_3)%*%OLh_3
  OLtheta <- matrix(0,2,164)
  OLF <- rep(list(matrix(0,2,2)),164)
  OLh <- matrix(0,2,164)
  for (j in (first):163){
    OLtheta[,first] <- OLtheta3
    OLF[[first]] <- OLF_3
    OLh[,first] <- OLh_3
    OLpred.nh[j+1] <-  c(1,1)%*%OLtheta[,j]
    OLF[[j+1]] <- OLF[[j]]+(k^j)*(c(1,-j))%*%t(c(1,-j))
    OLh[,j+1] <- (k)*LInv%*%OLh[,j]+(c(1,0))*Training_data$nh[j+1]
    OLtheta[,j+1] <- solve(OLF[[j+1]])%*%OLh[,j+1]
  }
  ## not using the first 3 values as they are not one-step predictions
  ## and not using the first 5 one-step predictions as they are used
  ## as the burn in period
  OLpred.err <- (Training_data$nh[9:164]-OLpred.nh[9:164])^2
  sum_sq_err[1000*k-844] <- sum(OLpred.err)
}
k <- c(seq(0.845,0.855,0.001))
plot(sum_sq_err~k,type="l",col="red",main="Sum of squared errors for values of lambda",
     xlab="Lambda",ylab="sum of squared errors")
## optimal lambda=0.848


## running the code from question 1.3 but with lambda=0.848 instead of
## lambda=0.8
first <- 3
IOLpred.nh <- numeric(164)
IOLpred.nh[1:3] <- NA
first_data <- c(Training_data$nh[1:3])
time <- Training_data$year-1849
L <- matrix(c(1,1,0,1),ncol=2)
LInv <- solve(L)

## creating the new intial epsilon matrix
IOlambda_init <- matrix(c((1/0.848^2),0,0, 0,(1/0.848),0, 0,0,1),nrow=3)


X3 <- cbind(1,(time[1:3])-first)
IOLF_3 <- t(X3)%*%solve(IOlambda_init)%*%X3
IOLh_3 <- t(X3)%*%solve(IOlambda_init)%*%first_data
IOLtheta3 <- solve(IOLF_3)%*%IOLh_3
IOLtheta <- matrix(0,2,164)
IOLF <- rep(list(matrix(0,2,2)),164)
IOLh <- matrix(0,2,164)

## one-step predicitions for the optimised local linear trend model
for (j in (first):163){
  IOLtheta[,first] <- IOLtheta3
  IOLF[[first]] <- IOLF_3
  IOLh[,first] <- IOLh_3
  IOLpred.nh[j+1] <-  c(1,1)%*%IOLtheta[,j]
  IOLF[[j+1]] <- IOLF[[j]]+(0.848^j)*(c(1,-j))%*%t(c(1,-j))
  IOLh[,j+1] <- (0.848)*LInv%*%IOLh[,j]+(c(1,0))*Training_data$nh[j+1]
  IOLtheta[,j+1] <- solve(IOLF[[j+1]])%*%IOLh[,j+1]
}

plot(IOLpred.nh~Training_data$year, type="l",
     main="Optimised local one-step predictions",
     xlab="Year",ylab="Temperature anomality",
     col="red",xlim=c(1850,2018),ylim=c(-0.668,1.064))
lines(Training_data$nh~Training_data$year,col="blue")
legend("topleft",legend=c("Training data",
                          "Optimised local one-step predictions"),
       col=c("blue","red"),lty=1:1, cex=0.8)

IOLpred.err <- Training_data$nh-IOLpred.nh  
IOLpred.err[1:3] <- NA
plot(IOLpred.err~Training_data$year,type="l",
     main="Optimised local one-step prediction errors",xlab="Year",
     ylab="Prediction error", col="red",xlim=c(1850,2018),
     ylim=c(-0.668,1.064))

IOLpred.test <- c(c(1,1)%*%IOLtheta[,164],c(1,2)%*%IOLtheta[,164],
                c(1,3)%*%IOLtheta[,164],c(1,4)%*%IOLtheta[,164],
                c(1,5)%*%IOLtheta[,164])

X <- cbind(1,-161:0)
IOLS <- (Training_data$nh[3:164]-X%*%IOLtheta[,164])
IOLepsilon <- matrix(0,162,162)
for (i in 162:1){
  IOLepsilon[163-i,163-i] <- 1/(0.848^(i-1))
}

## the divisor for IOLvar is T-p where T is 1/1-lambda, lambda=0.848
## and p=2
IOLvar <- (t(IOLS)%*%solve(IOLepsilon)%*%IOLS) / ((1/(1-0.848))-2)

IOLpredict.int <- numeric(5)
for(i in 1:5){
  IOLpredict.int[i] <-
    qt(0.975, 160)*sqrt(IOLvar)* 
    sqrt((1+t(c(1,i))%*%solve(IOLF[[164]])%*%c(1,i)))
}

IOLpred.test.lower <- IOLpred.test - IOLpredict.int
IOLpred.test.upper <- IOLpred.test + IOLpredict.int

plot(IOLpred.nh~Training_data$year, type="l",
     main="Optimised local one-step predictions with testing data predictions",
     xlab="Year",ylab="Temperature anomality",
     col="red",xlim=c(1850,2018),ylim=c(-0.668,1.064))
lines(Training_data$nh~Training_data$year,col="black")
lines(IOLpred.test~Testing_data$year,col="blue")
lines(IOLpred.test.lower~Testing_data$year,col="green")
lines(IOLpred.test.upper~Testing_data$year,col="green")
points(Testing_data$nh~Testing_data$year,pch=19,cex=0.5,col="orange")
legend("topleft",legend=c("Training data","Local one-step predictions",
                          "Testing data predictions",
                          "95% prediction interval","Testing data"),
       col=c("black","red","blue","green","orange"),lty=1:1, cex=0.8)



## table comparing the optimised local linear trend model 
## predictions of the testing data and the actual testing data.
##Also includes the upper and lower bounds from the prediction interval
IOLpredict.table <- data.frame(Testing_data$year,IOLpred.test,
                             IOLpred.test.lower,IOLpred.test.upper
                             ,Testing_data$nh)
names(IOLpredict.table) <- c("Year","Testing data predictions",
                           "95% prediction interval lower bound",
                           "95% prediction interval upper bound",
                           "Testing data")


## Question 1.5
## Quadratic trend model for lambda=0.848
first <- 4
QIOLpred.nh <- numeric(164)
QIOLpred.nh[1:4] <- NA
Qfirst_data <- c(Training_data$nh[1:3])
time <- Training_data$year-1849
QL <- matrix(c(1,1,1/2,0,1,1,0,0,1),ncol=3)
QLInv <- solve(QL)

## creating the new intial epsilon matrix
QIOlambda_init <- matrix(c((1/0.848^2),0,0, 0,(1/0.848),0, 0,0,1),nrow=3)


QX3 <- cbind(1,(time[1:3])-first,(((time[1:3])-first)^2)/2)
QIOLF_3 <- t(QX3)%*%solve(QIOlambda_init)%*%QX3
QIOLh_3 <- t(QX3)%*%solve(QIOlambda_init)%*%Qfirst_data
QIOLtheta3 <- solve(QIOLF_3)%*%QIOLh_3
QIOLtheta <- matrix(0,3,164)
QIOLF <- rep(list(matrix(0,3,3)),164)
QIOLh <- matrix(0,3,164)

## one-step predicitions for the optimised local quadratic trend model
for (j in (first):163){
  QIOLtheta[,first] <- QIOLtheta3
  QIOLF[[first]] <- QIOLF_3
  QIOLh[,first] <- QIOLh_3
  QIOLpred.nh[j+1] <-  c(1,1,1/2)%*%QIOLtheta[,j]
  QIOLF[[j+1]] <- QIOLF[[j]]+(0.848^j)*(c(1,-j,-(j^2)/2)%*%t(c(1,-j,-(j^2)/2)))
  QIOLh[,j+1] <- (0.848)*QLInv%*%QIOLh[,j]+(c(1,0,0))*Training_data$nh[j+1]
  QIOLtheta[,j+1] <- solve(QIOLF[[j+1]])%*%QIOLh[,j+1]
}

plot(QIOLpred.nh~Training_data$year, type="l",
     main="Optimised local one-step predictions for quadratic trend model",
     xlab="Year",ylab="Temperature anomality",
     col="red",xlim=c(1850,2018),ylim=c(-1,3.5))
lines(Training_data$nh~Training_data$year,col="blue")
legend("topright",legend=c("Training data",
                          "Optimised local one-step predictions
                          for quadratic trend model"),
       col=c("blue","red"),lty=1:1, cex=0.8)

QIOLpred.err <- Training_data$nh-QIOLpred.nh  
QIOLpred.err[1:4] <- NA
plot(QIOLpred.err~Training_data$year,type="l",
     main="Optimised local one-step prediction errors for quadratic trend model",
     xlab="Year",ylab="Prediction error", col="red",xlim=c(1850,2018),
     ylim=c(-3.5,1.064))

QIOLpred.test <- c(c(1,1,1/2)%*%QIOLtheta[,164],c(1,2,2)%*%QIOLtheta[,164],
                  c(1,3,9/2)%*%QIOLtheta[,164],c(1,4,8)%*%QIOLtheta[,164],
                  c(1,5,25/2)%*%QIOLtheta[,164])

v <- c(160:0)^2
QX <- cbind(1,-160:0,-v/2)
QIOLS <- (Training_data$nh[4:164]-QX%*%QIOLtheta[,164])
QIOLepsilon <- matrix(0,161,161)
for (i in 161:1){
  QIOLepsilon[162-i,162-i] <- 1/(0.848^(i-1))
}

## the divisor for QIOLvar is T-p where T is 1/1-lambda, lambda=0.848
## and p=3
QIOLvar <- (t(QIOLS)%*%solve(QIOLepsilon)%*%QIOLS) / ((1/(1-0.848))-3)

QIOLpredict.int <- numeric(5)
for(i in 1:5){
  QIOLpredict.int[i] <-
    qt(0.975, 160)*sqrt(QIOLvar)* 
    sqrt((1+t(c(1,i,(i^2)/2))%*%solve(QIOLF[[164]])%*%c(1,i,(i^2)/2)))
}

QIOLpred.test.lower <- QIOLpred.test - QIOLpredict.int
QIOLpred.test.upper <- QIOLpred.test + QIOLpredict.int

plot(QIOLpred.nh~Training_data$year, type="l",
     main="Optimised local one-step predictions with testing data predictions
     for quadratic trend model",
     xlab="Year",ylab="Temperature anomality",
     col="red",xlim=c(1850,2018),ylim=c(-4,3.5))
lines(Training_data$nh~Training_data$year,col="black")
lines(QIOLpred.test~Testing_data$year,col="blue")
lines(QIOLpred.test.lower~Testing_data$year,col="green")
lines(QIOLpred.test.upper~Testing_data$year,col="green")
points(Testing_data$nh~Testing_data$year,pch=19,cex=0.5,col="orange")
legend("bottomleft",legend=c("Training data","Local one-step predictions for quadratic trend model",
                          "Testing data predictions",
                          "95% prediction interval","Testing data"),
       col=c("black","red","blue","green","orange"),lty=1:1, cex=0.8)



## table comparing the optimised local quadratic trend model 
## predictions of the testing data and the actual testing data.
##Also includes the upper and lower bounds from the prediction interval
QIOLpredict.table <- data.frame(Testing_data$year,QIOLpred.test,
                               QIOLpred.test.lower,QIOLpred.test.upper
                               ,Testing_data$nh)
names(QIOLpredict.table) <- c("Year","Testing data predictions",
                             "95% prediction interval lower bound",
                             "95% prediction interval upper bound",
                             "Testing data")
