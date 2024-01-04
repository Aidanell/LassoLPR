library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)

lassoLPR <- function(x, y, bandwidth, evaluatePoints, numDerivatives = 10){
  
  lassoOutput <- matrix(nrow=401, ncol=numDerivatives+1) #Plus one for constant column
  
  for(i in 1:length(evaluatePoints)){
    point <- evaluatePoints[i]
    X <- buildFeatureMatrix(point, numTerms = numDerivatives, x)
    weights <- computeWeights(x, midpoint = point, bandwidth)
    
    lassoFit <- cv.glmnet(X, y, weights = weights, standardize=TRUE, alpha=1)
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    print(lassoCoef)
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  return(lassoOutput)
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
      else{weights[i] <- epan(diff)/ (bandwidth)}
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



plot(motorcycle)


motorcycle <- data.frame(read.csv("motorcycle.txt", head = TRUE, sep=" "))
colnames(motorcycleData)

bandwidth <- dpill(motorcycle[,1], motorcycle[,2])
lassoBandwidth <- bandwidth * 8
pointstoEvaluateAt <- seq(min(motorcycle[,1]), max(motorcycle[,1]), length.out=401)


motorLasso <- lassoLPR(motorcycle[,1], motorcycle[,2], lassoBandwidth, pointstoEvaluateAt)
lines(pointstoEvaluateAt, motorLasso[,1], col='red')

motorLocpoly <- locpoly(motorcycle[,1], motorcycle[,2], bandwidth = bandwidth, degree=1)
lines(motorLocpoly)



###First Derivative
motorLocPoly1st <- locpoly(motorcycle[,1], motorcycle[,2], bandwidth = bandwidth, drv=1)
plot(motorLocPoly1st, type='l')
lines(pointstoEvaluateAt, motorLasso[,2], col='red')






