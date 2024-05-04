library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)

#Original Method----
lassoLPR <- function(x, y, bandwidth, evaluatePoints, numDerivatives = 10){
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
      else{weights[i] <- dnorm(diff)/bandwidth}
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

#Testing Method with lambda 0----
Testinglambda0 <- function(x, y, bandwidth, evaluatePoints, numDerivatives = 10){
  
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
      else{weights[i] <- dnorm(diff)/bandwidth}
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
  
  
  lassoOutput <- matrix(nrow=401, ncol=numDerivatives+1) #Plus one for constant column
  ListOfLambdas <- c()
  for(i in 1:length(evaluatePoints)){
    point <- evaluatePoints[i]
    X <- buildFeatureMatrix(point, numTerms = numDerivatives, x)
    weights <- computeWeights(x, midpoint = point, bandwidth)
    
    lassoFit <- glmnet(X, y, weights = weights, alpha=1)
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s=0)
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
    
  }
  
  lassoData <- vector(mode='list', 2)
  lassoData[[1]] <- lassoOutput
  return(lassoData)
}

#Testing getting rid of factorial----
Testing <- function(x, y, bandwidth, evaluatePoints, numDerivatives = 10){
  
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
      else{weights[i] <- dnorm(diff)/bandwidth}
    }
    
    return(weights)
  }
  
  
  
  buildFeatureMatrix <- function(midpoint, numTerms, xValues){
    # Every j'th column is the coefficient of the (j-1)'th derivative in the
    # Taylor polynomial formed at the i'th point.
    
    X <- matrix(nrow = length(xValues), ncol = numTerms)
    
    for(i in 1:numTerms){
      X[,i] <- ((xValues - midpoint)^(i))
    }
    return(X)
  }
  
  
  lassoOutput <- matrix(nrow=401, ncol=numDerivatives+1) #Plus one for constant column
  ListOfLambdas <- c()
  for(i in 1:length(evaluatePoints)){
    point <- evaluatePoints[i]
    X <- buildFeatureMatrix(point, numTerms = numDerivatives, x)
    weights <- computeWeights(x, midpoint = point, bandwidth)
    
    lassoFit <- glmnet(X, y, weights = weights, alpha=1)
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s=0)
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef)
    
    for(i in 1:401){
      for(j in 1:numDerivatives+1){
        lassoOutput[i,j] <- lassoOutput[i,j] / factorial(j-1)
      }
    }
    
  }
  
  lassoData <- vector(mode='list', 2)
  lassoData[[1]] <- lassoOutput
  return(lassoData)
}






#Evaluating how lambas change output----
lambdaTesting <- function(x, y, bandwidth, evaluatePoints, numDerivatives = 10){
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
      else{weights[i] <- dnorm(diff)/bandwidth}
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
  
  lassoOutput <- matrix(nrow=401, ncol=numDerivatives+1) #Plus one for constant column
  
    X <- buildFeatureMatrix(evaluatePoints, numTerms = numDerivatives, x)
    weights <- computeWeights(x, midpoint = evaluatePoints, bandwidth)
    
    lassoFit <- glmnet(X, y, weights = weights, standardize=TRUE, alpha=1)
    
    return(list(lassoFit, X))
} 
##Script----
#Getting Data
motorcycleData <- data.frame(read.csv("motorcycle.txt", head = TRUE, sep=" "))
colnames(motorcycleData)

#Bandwidth and evalPoints
bandwidth <- dpill(motorcycle[,1], motorcycle[,2])
lassoBandwidth <- bandwidth * 8
pointstoEvaluateAt <- seq(min(motorcycle[,1]), max(motorcycle[,1]), length.out=401)


#Evaluate lasso/locpoly methods
motorLasso <- lassoLPR(motorcycle[,1], motorcycle[,2], lassoBandwidth, pointstoEvaluateAt, 10)
motorLocpoly <- locpoly(motorcycle[,1], motorcycle[,2], bandwidth = bandwidth, degree=1)
#motorTesting <- Testing(motorcycle[,1], motorcycle[,2], lassoBandwidth, pointstoEvaluateAt)

#Plotting
par(mfrow=c(1,1))
plot(motorcycleData, main='Lasso vs Locpoly')
lines(motorLocpoly)
lines(pointstoEvaluateAt, motorLasso[[1]][,1], col='red')
#lines(pointstoEvaluateAt, motorTesting[[1]][,1], col='darkgreen')
# 
#plot(pointstoEvaluateAt, motorLasso[[2]], pch=20)
