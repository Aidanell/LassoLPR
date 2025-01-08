
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
p <- 4
sigma <- sqrt(0.5)
deriv <- 3 #Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,gaussK)


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

x <- sampleData$x
trueY <- eval(peak)

testingCustomCV <- function(x, y, h, p=10, numGridPoints=401){
  
  gridPoints <- seq(min(x), max(x), length.out=numGridPoints)
  
  
  lassoOutput <- matrix(nrow=numGridPoints, ncol=p+1) #ith column shows estimates of (i-1)th derivative 
  lambdas <- numeric(numGridPoints) #Tracks lambda parameter over x
  
  
  #Estimate lasso polynomial at every gridpoint
  for(i in 1:numGridPoints){
    
    currentPoint <- gridPoints[i]
    X <- buildFeature(currentPoint, p, x)
    lassoWeights <- computeWeights(x, currentPoint, h, kernel='norm')
    hist(lassoWeights)
    closeIndex <- which(lassoWeights > quantile(lassoWeights, 0.5)) #I do not care about performance if they are not in the window
    #closeIndex <- which( abs(x - currentPoint)/1.3*h < 1)
    print(length(closeIndex))
    
    lambdaSeq <- glmnet::glmnet(X[closeIndex,], y[closeIndex], weights = lassoWeights[closeIndex], maxit=10**7, nlambda=100)$lambda

    PerformCV <- NfoldCV(x[closeIndex], y[closeIndex], currentPoint, lambdaSeq, lassoWeights[closeIndex], p, trueY[closeIndex], folds=20)
    lassoCoef <- PerformCV[[1]]
    lambdas[i] <- PerformCV[[2]]
    print(i)
    #The i'th row contains all estimated p+1 derivatives estimation at the i'th gridpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  return(list(lassoOutput, lambdas))
}

testLasso <- testingCustomCV(sampleData$x, sampleData$y, h=bandwidth, p=9)


#Choose what derivative to plot, Change as needed*********
derivative <- 1

#Runs a locpoly example to compare with.
locBand <- dpill(sampleData$x, sampleData$y)
locOut <- locpoly(sampleData$x, sampleData$y, bandwidth=locBand*2, drv=derivative)

#Plotting
title <- paste("Performance Comparison of Lasso vs Locpoly for Derivative ", derivative+1)
plot(seq(0,1, length.out=401), trueFunctions[[derivative+1]], type='l', lwd=2, ylab='y', xlab='x', main=title)
lines(locOut, col='red', lwd=2)
lines(seq(0,1, length.out=401), testLasso[[1]][,derivative+1], col='blue', type='l', lwd=2)
#Legend
legend('bottomleft', legend=c('LassoLPR', "Locpoly", "True Function"), fill=c("blue", "red", "black"))
#Optionally plot with Datapoints (only relevant for derivative=0)
if(derivative == 0){
  points(sampleData$x, sampleData$y, pch=20, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1))  
}

plot(seq(0,1, length.out=401), testLasso[[1]][,2], type='l')
plot(seq(0,1, length.out=401), testLasso[[2]], type='l')


