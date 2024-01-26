library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)

lassoLPR <- function(x, y, bandwidth, evaluatePoints, numDerivatives = 10){
  
  lassoOutput <- matrix(nrow=401, ncol=numDerivatives+1) #Plus one for constant column
  ListOfLambdas <- c()
  for(i in 1:length(evaluatePoints)){
    point <- evaluatePoints[i]
    X <- buildFeatureMatrix(point, numTerms = numDerivatives, x)
    weights <- computeWeights(x, midpoint = point, bandwidth)
    
    lassoFit <- cv.glmnet(X, y, weights = weights, standardize=FALSE, alpha=1)
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





motorcycleData <- data.frame(read.csv("motorcycle.txt", head = TRUE, sep=" "))
colnames(motorcycleData)

bandwidth <- dpill(motorcycle[,1], motorcycle[,2])
lassoBandwidth <- bandwidth * 4
lassoBandwidth
pointstoEvaluateAt <- seq(min(motorcycle[,1]), max(motorcycle[,1]), length.out=401)


motorLasso <- lassoLPR(motorcycle[,1], motorcycle[,2], lassoBandwidth, pointstoEvaluateAt)
#lines(pointstoEvaluateAt, motorLasso[,1], col='red')

motorLocpoly <- locpoly(motorcycle[,1], motorcycle[,2], bandwidth = bandwidth, degree=1)
#plot(motorcycleData)
#lines(motorLocpoly)


par(mfrow=c(1,2))
plot(seq(min(motorcycle[,1]), max(motorcycle[,1]), length.out=401), motorLasso[[2]], main='Lambda Magnitude Over Time')
plot(motorcycleData, main='Lasso vs Locpoly')
lines(motorLocpoly)
lines(pointstoEvaluateAt, motorLasso[[1]][,1], col='red')

# motorLasso[[2]]
# ###First Derivative
# motorLocPoly1st <- locpoly(motorcycle[,1], motorcycle[,2], bandwidth = bandwidth, drv=1)
# plot(motorLocPoly1st, type='l')
# lines(pointstoEvaluateAt, motorLasso[[1]][,2], col='red')
# 
# 
# motorLocPoly2nd <- locpoly(motorcycle[,1], motorcycle[,2], bandwidth = bandwidth, drv=2)
# plot(pointstoEvaluateAt, motorLasso[[1]][,3], col='red', type='l')
# lines(motorLocPoly2nd, type='l')
#
motorLasso[[1]]
