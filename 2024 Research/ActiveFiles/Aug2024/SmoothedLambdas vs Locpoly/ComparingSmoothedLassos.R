library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)
library(locpol)


left <- 0
right <- 1

numTerms <- 10
numPoints = 100
set.seed(100)
bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2)) #sigma = 0.1
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2)) # sigma = sqrt(0.5)
sine <- expression(sin(5*pi*x)) #sigma = 0.5
evalPoints <- seq(0,1, length.out=401)

equation <- peak
sigma <- sqrt(0.5)

#Calculates all the derivatives needed for plotting exact functions and calculating errors
derivList <- derivCalc(func = equation, numDeriv = numTerms)

##############Building Simulated Data##################
x <- sort(runif(numPoints, min=left, max=right))
noise <- rnorm(n = length(x), mean = 0, sd = sigma)
exactY <- eval(equation)
y <- exactY + noise





#Test 1
#Bandwidth = Dpill * 8, Epan Weights
set.seed(1)
bandwidth1 <- dpill(x,y)*8
lassoTest1 <- LassoSmoothed(x,y, bandwidth1, TRUE, numTerms)


#Test 2
#Bamdwidth = Thumbbw(degree=0, epaK), Epan Weights
set.seed(1)
bandwidth2 <- thumbBw(x,y,0,EpaK)
lassoTest2 <- LassoSmoothed(x,y,bandwidth2, TRUE, numTerms)

#Test 3
#Bandwidth = Thumbbw(degree=0, gaussK), Normal Weights
set.seed(1)
bandwidth3 <- thumbBw(x,y,0,gaussK)
lassoTest3 <- LassoSmoothed(x,y,bandwidth3, FALSE, numTerms)

#Test 4
#Bamdwidth = Thumbbw(degree=3, epaK), Epan Weights
set.seed(1)
bandwidth4 <- thumbBw(x,y,3,EpaK)
lassoTest4 <- LassoSmoothed(x,y,bandwidth4, TRUE, numTerms)

#Test 5
#Bamdwidth = Thumbbw(degree=1, epaK), Epan Weights
set.seed(1)
bandwidth5 <- thumbBw(x,y,1,EpaK)
lassoTest5 <- LassoSmoothed(x,y,bandwidth5, TRUE, numTerms)

#Test 6
#Bamdwidth = Thumbbw(degree=3, epaK), Epan Weights
set.seed(1)
bandwidth6 <- thumbBw(x,y,3,EpaK)
lassoTest6 <- LassoSmoothed(x,y,bandwidth6, TRUE, numTerms)

###Plotting
temp <- x
x <- evalPoints
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),eval(derivList[6]))
x <- temp


#Plotting Graphs
deriv <- 0
plot(evalPoints, trueFunctions[[deriv+1]], type='l', lwd=2)
lines(evalPoints, lassoTest1[[2]][,deriv+1], col='blue', lwd=2)
#lines(evalPoints, lassoTest2[[2]][,deriv+1], col='orange', lwd=2)
#lines(evalPoints, lassoTest5[[2]][,deriv+1], col='purple', lwd=2)
lines(evalPoints, lassoTest6[[2]][,deriv+1], col="green", lwd=2)
#lines(evalPoints, lassoTest2[[2]][,deriv+1], col='red', lwd=2)
lines(evalPoints, lassoTest3[[2]][,deriv+1], col='red', lwd=2)


#Plotting  Smooth Lambdas
plot(evalPoints,lassoTest2[[4]], type='l', col='orange')
lines(evalPoints,lassoTest5[[4]], type='l', col='purple')
lines(evalPoints,lassoTest6[[4]], type='l', col="green")

#Plotting UnSmooth Lambdas
plot(evalPoints,lassoTest2[[3]], col='orange')
points(evalPoints,lassoTest5[[3]], col='purple')
points(evalPoints,lassoTest6[[3]], col="green")


#Locpoly Calculations
deriv <- 4
locpolyBand <- thumbBw(x,y,deriv+1, gaussK)
locpolyFit <- KernSmooth::locpoly(x,y,degree=10,drv=deriv, bandwidth=locpolyBand, gridsize=401)


#Plotting comparing with Locpoly
title <- paste(c("Derivative", deriv, "p=10, n=100"))
plot(evalPoints, trueFunctions[[deriv+1]], type='l', lwd=2, main=title)
lines(evalPoints, lassoTest6[[2]][,deriv+1], col="orange", lwd=2)
lines(locpolyFit, col='purple', lwd=2)
legend("bottomleft", legend=c("True Function", "SmoothedLasso", "Locpoly"), col=c('black',"orange","purple"),lty=1, lwd=2)




#LOOCV comparison
deriv <- 5
title <- paste(c("Derivative", deriv, "p=10, n=1000"))
plot(evalPoints, trueFunctions[[deriv+1]], type='l', lwd=2, main=title)
lines(evalPoints, lassoTest6[[1]][,deriv+1], col="orange", lwd=2)
lines(evalPoints, CV1000Lasso[[1]][,deriv+1], col="lightgreen", lwd=2)
