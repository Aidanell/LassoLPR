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
    
    lassoFit <- cv.glmnet(X, y, weights = weights, standardize=TRUE, alpha=1, maxit=10**7)
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
  MADEData0 <- matrix(0, 401, 1)
  MADEData1 <- matrix(0, 401, 1)
  
  #Script parameters
  n <- 200
  sigma <- 0.5
  equation <- expression(x + 2 * exp(-16*(x**2)))
  derivative <- expression(1 - 64*x*exp(-16*(x**2)))
  left <- -2
  right <- 2
  
  
  #Getting x and y values of gridpoints set
  x <- seq(left, right, length.out=401)
  gridpointsY <- eval(equation)
  gridDerivativeY <- eval(derivative)
  gridpointsX <- x
  
  
  for(i in 1:numIter){
    #Preparing simulated data                     
    x <- sort(runif(n, min=left, max=right))
    noise <- rnorm(n = length(x), mean = 0, sd = sigma)
    exactY <- eval(equation)
    y <- exactY + noise
    
    #Bandwidths
    bandwidth <- dpill(x,y)
    lassoBandwidth <- bandwidth * 8
    
    #Computing Non-parametric regression techniques
    lassoOutput <- lassoLPR(x, y, lassoBandwidth, gridpointsX)
    locpolyOutput <- locpoly(x, y, bandwidth = bandwidth, degree=1)
    
    #Computing Ratio data
    LassoMADE0 <- abs(lassoOutput[[1]][,1] - gridpointsY)
    LocpolyMADE0 <- abs(as.numeric(unlist(locpolyOutput$y)) - gridpointsY)
    MADEratio0 <- LocpolyMADE0 / LassoMADE0
    
    
    
    locpolyDer <- locpoly(x, y, bandwidth = bandwidth, degree=1, drv=1)
    
    LassoMADE1 <- abs(lassoOutput[[1]][,2] - gridDerivativeY)
    LocpolyMADE1 <- abs(as.numeric(unlist(locpolyDer$y)) - gridDerivativeY)
    MADEratio1 <- LocpolyMADE1 / LassoMADE1
    
    
    #Plotting
    par(mfrow=c(2,2))
    
    plot(x,exactY, type='l', main="Regression")
    lines(gridpointsX, lassoOutput[[1]][,1], col='red')
    lines(locpolyOutput, col='blue')
    points(x,y, pch=20)
    plot(gridpointsX, MADEratio0, ylim=c(0,5), main='Ratio of MADE', pch=20)
    abline(h=1, col='darkgreen')
    
    curve(1 - 64*x*exp(-16*(x**2)), from=-2, to=2, main="Regression of First Derivative")
    lines(gridpointsX, lassoOutput[[1]][,2], col='red')
    lines(locpolyDer, col='blue')
    plot(gridpointsX, MADEratio1, ylim=c(0,5), main="Ratio of MADE for First Derivative", pch=20)
    abline(h=1, col='darkgreen')
    
    MADEData0 <- cbind(MADEData0, MADEratio0)
    MADEData1 <- cbind(MADEData1, MADEratio1)

    MADEData <- vector(mode="list", length=2)
    MADEData[[1]] <- MADEData0
    MADEData[[2]] <- MADEData1
    
    print(i)
  }
  return(MADEData)
}

MADESim <- runMADESimulations(2)



gridpoints <- seq(-2, 2, length.out=401)

par(mfrow=c(2,2))
plot(gridpoints, apply(MADESim[[1]], 1, median, na.rm=TRUE), pch=20, main="Median MADE ratio")
abline(h=1, col='darkgreen')
plot(gridpoints, apply(MADESim[[2]], 1, median, na.rm=TRUE), pch=20, main="Median MADE ratio of First Derivative")
abline(h=1, col='darkgreen')

plot(gridpoints, apply(MADESim[[1]], 1, mean, na.rm=TRUE), pch=20, main="Average MADE ratio", ylim=c(0,10))
abline(h=1, col='darkgreen')
plot(gridpoints, apply(MADESim[[2]], 1, mean, na.rm=TRUE), pch=20, main="Average MADE ratio of First Derivative", ylim=c(0,10))
abline(h=1, col='darkgreen')







