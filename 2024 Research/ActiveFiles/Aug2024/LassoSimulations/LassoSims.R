####Creating a large scale Monte-Carlo Simulation
library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)
library(locpol)

#Parameters
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2)) # sigma = sqrt(0.5)

epochs <- 10

equation <- peak
n <- 500
p <- 10
sigma <- sqrt(0.5)
gridPoints <- seq(0,1,length.out=401)
deriv <- 1 #Used for bandwidth calculations

MasterDataList <- vector("list", epochs)
for(i in 1:epochs){
  MasterDataList[[i]] <- lassoEpoch(n, p, equation, sigma, deriv)
  print(i)
}

#True Output
derivList <- derivCalc(func = equation, numDeriv = 5)
x <- gridPoints
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),eval(derivList[6]))


#Plotting Average Output
aveOutput <- calcAverage(MasterDataList)
plot(gridPoints, trueFunctions[[1]])
lines(gridPoints, aveOutput$aveUnsmoothOutput, col='red')
lines(gridPoints, aveOutput$aveSmoothOutput, col='orange')
lines(gridPoints, aveLocOut, col='blue')

#Plotting Bias
plot(gridPoints, abs(aveOutput$aveUnsmoothOutput - trueFunctions[[1]]), col='red', type='l', lwd=2)
lines(gridPoints, abs(aveOutput$aveSmoothOutput- trueFunctions[[1]]), col='orange', lwd=2)
lines(gridPoints, abs(aveLocOut - trueFunctions[[1]]), col='blue', lwd=2)

#Plotting Variance
varOutput <- calcVariance(MasterDataList)
plot(gridPoints, varOutput$varUnsmoothOutput, col='red', type='l', lwd=2)
lines(gridPoints, varOutput$varSmoothOutput, col='orange', type='l', lwd=2)


calcVariance <- function(MasterData, deriv=0){
  epochs <- length(MasterData)
  
  varUnsmoothOutput <- rep(0, 401)
  varSmoothOutput <- rep(0, 401)
  
  aveOutput <- calcAverage(MasterData, deriv=deriv)
  for(i in 1:epochs){
    varUnsmoothOutput <- rbind(varUnsmoothOutput, (MasterData[[i]]$lasso[,deriv+1] - aveOutput$aveUnsmoothOutput)**2)
    varSmoothOutput <- rbind(varSmoothOutput, (MasterData[[i]]$smoothLasso[,deriv+1] - aveOutput$aveSmoothOutput)**2)
  }
  return(list(varUnsmoothOutput = colMeans(varUnsmoothOutput), varSmoothOutput=colMeans(varSmoothOutput)))
  
}

calcAverage <- function(MasterData, deriv = 0){
  epochs <- length(MasterData)
  
  aveUnsmoothOutput <- rep(0, 401)
  aveSmoothOutput <- rep(0, 401)
  for(i in 1:epochs){
    aveUnsmoothOutput <- rbind(aveUnsmoothOutput, MasterData[[i]]$lasso[,deriv+1])
    aveSmoothOutput <- rbind(aveSmoothOutput, MasterData[[i]]$smoothLasso[,deriv+1])
  }
  return(list(aveUnsmoothOutput = colMeans(aveUnsmoothOutput), aveSmoothOutput=colMeans(aveSmoothOutput)))
}

lassoEpoch <- function(n, p, equation, sigma, deriv){
  
  data <- buildData(n, equation, sigma)
  bandwidth <- thumbBw(data$x,data$y,deriv,EpaK)
  lassoData <- lassoCalculations(data$x, data$y, p, bandwidth)
  return(lassoData)
}


lassoCalculations <- function(x,y,p,bandwidth){
  
  gridPoints <- seq(0, 1, length.out=401)
  lassoOutput <- matrix(nrow=401, ncol=p)
  ListOfLambdas <- c()
  
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    X <- buildFeature(currentPoint, p, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    lassoFit <- cv.glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
    ListOfLambdas <- c(ListOfLambdas, lassoFit$lambda.min)
    
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  
  
  smoothedLambas <- lowess(gridPoints, ListOfLambdas, f=1/10)$y
  smoothLassoOutput <- matrix(nrow=401, ncol=p)
  
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    
    X <- buildFeature(currentPoint, p, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    lassoFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
    
    # Selects the desired coefficients with penalty lambda
    smoothLassoCoef <- coef(lassoFit, s=smoothedLambas[i])
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef) 
  }
  
  return(list(lasso=lassoOutput, smoothLasso=smoothLassoOutput,
              unsmoothedLambdas=ListOfLambdas, smoothedLambdas=smoothedLambas))
}


buildData <- function(n, equation, sigma){
  x <- sort(runif(n, min=0, max=1))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  trueY <- eval(equation)
  y <- trueY + noise
  return(list(x=x,y=y,trueY=trueY))
}


#Helper Methods----
# This function builds our feature matrix.
buildFeature <- function(midpoint, p, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = p-1)
  
  for(i in 1:p-1){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}


#Calculates true derivatives needed to compare estimates of derivatives with real
derivCalc <- function(func, numDeriv){
  drvList <- c(func)
  nextDrv <- func
  for(i in 1:numDeriv){
    nextDrv <- D(nextDrv, 'x')
    drvList <- c(drvList, nextDrv)
  }
  return(drvList)
}


#Function of the Epanechnikov kernel
epan <- function(x, bandwidth){
  return(3/4 * (1-x**2))
}


#Computes the scaled kernel weights for every input x
computeWeights <- function(x, midpoint, bandwidth, epan=TRUE){
  output <- vector(length=length(x))
  
  for(i in 1:length(x)){
    
    diff <- abs(x[i]-midpoint) / bandwidth
    
    if(epan){
      if(diff > 1){output[i] <- 0}
      else{output[i] <- epan(diff)/ (bandwidth)}
    }else{
      output[i] <- dnorm(diff)
    }
    
  }
  return(output)
}

