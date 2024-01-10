library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)

#Helper Methods----
lassoLPR <- function(x, y, bandwidth, evaluatePoints, numDerivatives = 10){
  
  lassoOutput <- matrix(nrow=401, ncol=numDerivatives+1) #Plus one for constant column
  ListOfLambdas <- c()
  for(i in 1:length(evaluatePoints)){
    point <- evaluatePoints[i]
    X <- buildFeatureMatrix(point, numTerms = numDerivatives, x)
    weights <- computeWeights(x, midpoint = point, bandwidth)
    
    lassoFit <- cv.glmnet(X, y, weights = weights, standardize=TRUE, alpha=1)
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    ListOfLambdas <- c(ListOfLambdas, lassoFit$lambda.min)
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
    
  }
  
  lassoData <- vector(mode='list', 2)
  lassoData[[1]] <- lassoOutput
  lassoData[[2]] <- ListOfLambdas
  return(lassoData)
} 


#Computes the scaled kernel weights for every input x
computeWeights <- function(x, midpoint, bandwidth){
  
  #Function of the Epanechnikov kernel
  epan <- function(x, bandwidth){
    return(3/4 * (1-x**2))
  }
  
  weights <- vector(length=length(x))
  for(i in 1:length(x)){
    
    diff <- abs(x[i]-midpoint) / bandwidth
    if(diff > 1){weights[i] <- 0}
    else{weights[i] <- dnorm(diff)/ (bandwidth)}
  }
  
  return(weights)
}



buildFeatureMatrix <- function(midpoint, numTerms, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = numTerms)
  
  for(i in 1:numTerms){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}

runMADESimulations <- function(numIter){
  MADEData <- matrix(0, 401, 1)
  
  #Script parameters
  n <- 200
  sigma <- 0.5
  equation <- expression(x + 2 * exp(-16*(x**2)))
  left <- -2
  right <- 2
  
  
  #Getting x and y values of gridpoints set
  x <- seq(left, right, length.out=401)
  gridpointsY <- eval(equation)
  gridpointsX <- x
  
  #Bandwidths
  bandwidth <- dpill(x,y)
  lassoBandwidth <- bandwidth * 8
  
  for(i in 1:numIter){
    #Preparing simulated data                     
    x <- sort(runif(numPoints, min=left, max=right))
    noise <- rnorm(n = length(x), mean = 0, sd = sigma)
    exactY <- eval(equation)
    y <- exactY + noise
    
    #Computing Non-parametric regression techniques
    lassoOutput <- lassoLPR(x, y, lassoBandwidth, gridpointsX)
    locpolyOutput <- locpoly(x, y, bandwidth = bandwidth, degree=1)
    
    #Computing Ratio data
    LassoMADE <- abs(lassoOutput[[1]][,1] - gridpointsY)
    LocpolyMADE <- abs(as.numeric(unlist(locpolyOutput$y)) - gridpointsY)
    MADEratio <- LassoMADE / LocpolyMADE
    
    
    #Plotting
    par(mfrow=c(1,2))
    plot(x,exactY, type='l')
    lines(gridpoints, lassoOutput[[1]][,1], col='red')
    lines(locpolyOutput, col='blue')
    plot(gridpointsY, MADEratio, ylim=c(0,5))
    
    
    MADEData <- cbind(MADEData, MADEratio)
    print(nrow(MADEData))
    print(length(MADEratio))
  }
  return(MADEData)
}

MADESim <- runMADESimulations(2)

MADESim
par(mfrow=c(1,1))
plot(gridpointsY, apply(MADESim, 1, median, na.rm=TRUE),type='l', ylim=c(0,2))

sum(apply(MADESim, 1, median, na.rm=TRUE) < 1)

apply(MADESim, 1, median, na.rm=TRUE)
