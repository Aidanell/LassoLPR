
library(locpol)
library(KernSmooth)
#This file is to test to try and find a hypothetical optimal lambda GIVEN we know the true answer
#If one is not found that performs ideal for a large amount of derivatives, it may be an issue
set.seed(20)

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
deriv <- 0 #Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,gaussK)
gridPoints <- seq(0, 1, length.out=401)

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


deriv <- 0
lamOut <- numeric(length(gridPoints))
predOut <- numeric(length(gridPoints))

findOptimalLambda <- function(deriv, i){
  
  
  #Find good sequence of lambdas
  #Fit for all lambdas and choose one which optimizes against the first p derivatives
  #Weight based on what derivative we want?
  
  X <- buildFeature(gridPoints[i], p, sampleData$x) 
  lassoWeights <- computeWeights(sampleData$x, gridPoints[i], bandwidth, kernel='norm') 
  
  
  lassoFit <- glmnet::glmnet(X, sampleData$y, weights = lassoWeights, maxit=10**7, nlambda=50)
  lambdas <- lassoFit$lambda
  
  MAE <- numeric(length(lambdas))
  predicted <- numeric(length(lambdas))
  
  truePoint <- trueFunctions[[deriv+1]][i]
  
  for(j in 1:length(lambdas)){
    predicted[j] <- as.vector(coef(lassoFit, s=lambdas[j]))[deriv+1]
    MAE[j] <- abs(predicted[j] - truePoint)
  }
  
  
  #Calculating minimum
  #Note: If no error at multiple lambdas, which.min would choose the larger lambda. We do not want that
  
  possibleMin <- which(MAE == min(MAE))
  lastMinIndex <- tail(possibleMin, 1)
  minLam <- lambdas[lastMinIndex]
  plot(lambdas, MAE)
  points(minLam, MAE[lastMinIndex], pch=16, cex=2, col='red')
  
  return(c(minLam, predicted))
}

for(i in 1:length(gridPoints)){
  
  currOut <- findOptimalLambda(deriv, i)
  lamOut[i] <- currOut[1]
  predOut[i] <- currOut[2]
  print(i)
}



locBand <- dpill(sampleData$x, sampleData$y)
locOut <- locpoly(sampleData$x, sampleData$y, bandwidth=locBand*(deriv+1), drv=deriv)

plot(gridPoints, lamOut, type='l')


plot(gridPoints, trueFunctions[[deriv+1]], lwd=2, type='l')
lines(gridPoints, predOut, type='l', col='blue', lwd=2)
lines(locOut, col='red', lwd=2)

#Why can't lambda get all the way to 0? So no chosen lambda is big enough to bring it down
#However, glmnet chooses lambda based on the value that would completely zero out everything
#Therefore, I think the lambda is too large 

#So if we want to estimate a certain lambda, we are always going to be restricted
#by how glmnet chooses variables


