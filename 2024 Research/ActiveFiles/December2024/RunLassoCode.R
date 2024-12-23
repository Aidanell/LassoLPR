library(locpol)
library(KernSmooth)


#----Pre-Setup----
#Function that builds simulated data
buildData <- function(n, equation, sigma){
  x <- sort(runif(n, min=0, max=1))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  trueY <- eval(equation)
  y <- trueY + noise
  return(list(x=x,y=y,trueY=trueY))
}

#Create neccesary data----
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5*x)
p <- 9
sigma <- sqrt(0.5)
deriv <- 1 #Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,EpaK)


#True Functions----
derivCalc <- function(func, numDeriv){
  drvList <- c(func)
  nextDrv <- func
  for(i in 1:numDeriv){
    nextDrv <- D(nextDrv, 'x')
    drvList <- c(drvList, nextDrv)
  }
  return(drvList)
}

#Builds the true datapoints for the function. Helps with MSE
derivList <- derivCalc(func = peak, numDeriv = 5)
x <- seq(0,1,length.out=401)
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),eval(derivList[6]))


##--Actually Performing Test----

#This runs LassoLPR from the other file.
testLasso <- lassoLPR(sampleData$x, sampleData$y, h=bandwidth, p=p)

#Choose what derivative to plot, Change as needed*********
derivative <- 0

#Runs a locpoly example to compare with.
locBand <- dpill(sampleData$x, sampleData$y)
locOut <- locpoly(sampleData$x, sampleData$y, bandwidth=locBand, drv=derivative)

#Plotting
title <- paste("Performance Comparison of Lasso vs Locpoly for Derivative ", derivative+1)
plot(seq(0,1, length.out=401), testLasso$lasso[,derivative+1], type='l', col='blue', lwd=2, ylab='y', xlab='x', main=title)
lines(locOut, col='red', lwd=2)
lines(x, trueFunctions[[derivative+1]], type='l', lwd=2)
#Legend
legend('bottomleft', legend=c('LassoLPR', "Locpoly", "True Function"), fill=c("blue", "red", "black"))

#Optionally plot with Datapoints (only relevant for derivative=0)
if(derivative == 0){
  points(sampleData$x, sampleData$y, pch=20, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1))  
}


