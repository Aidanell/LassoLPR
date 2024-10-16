testLambdasLassoLPR <- function(x, y, p, bandwidth){
  
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
  
  #Method Body----
  gridPoints <- seq(0, 1, length.out=401)
  lassoOutput <- matrix(nrow=401, ncol=p)
  ListOfLambdas <- c()
  
  #Estimate lasso polynomial at every gridpoint
  for(i in 1:length(gridPoints)){
    
    currentPoint <- gridPoints[i]
    X <- buildFeature(currentPoint, p, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    
    #Which lambdas to try?
    if(i > 1){
        lastLambda <- ListOfLambdas[i-1]
        numLambda <- 15
        testLambdas <- seq(lastLambda-0.005, lastLambda+0.005, length.out=numLambda)
        testLambdas <- ifelse(testLambdas < 0, 0, testLambdas)
        print(testLambdas)
        lassoFit <- cv.glmnet(X, y, weights = lassoWeights, maxit=10**7, lambda=testLambdas)
        
    }else{
      lassoFit <- cv.glmnet(X, y, weights = lassoWeights, maxit=10**7)
      }
    ListOfLambdas <- c(ListOfLambdas, lassoFit$lambda.min)
    
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  
  #Smooth Lambdas output with loess, and refit with new lambda graph. This reduces variance
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
  
  #Return list of smooth and unsmooth lasso, alongside the lambdas
  return(list(lasso=lassoOutput, smoothLasso=smoothLassoOutput,
              unsmoothedLambdas=ListOfLambdas, smoothedLambdas=smoothedLambas))
}


averageLassoLPR <- function(x, y, p, bandwidth){
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
        output[i] <- dnorm(diff)/bandwidth
      }
      
    }
    return(output)
  }
  
  #Method Body----
  gridPoints <- seq(0, 1, length.out=401)
  lassoOutput <- matrix(nrow=401, ncol=p)
  ListOfLambdas <- c()
  
  #Estimate lasso polynomial at every gridpoint
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
  
  #Smooth Lambdas output with loess, and refit with new lambda graph. This reduces variance
  smoothedLambas <- rep(mean(ListOfLambdas), 401)
  smoothLassoOutput <- matrix(nrow=401, ncol=p)
  
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    
    X <- buildFeature(currentPoint, p, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    lassoFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
    
    # Selects the desired coefficients with penalty lambda
    smoothLassoCoef <- coef(lassoFit, s=mean(smoothedLambas))
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef) 
  }
  
  #Return list of smooth and unsmooth lasso, alongside the lambdas
  return(list(lasso=lassoOutput, smoothLasso=smoothLassoOutput,
              unsmoothedLambdas=ListOfLambdas, smoothedLambdas=smoothedLambas))
}


excludeZeroWeightsLassoLPR <- function(x, y, p, bandwidth){
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
        output[i] <- dnorm(diff)/bandwidth
      }
      
    }
    return(output)
  }
  
  #Method Body----
  gridPoints <- seq(0, 1, length.out=401)
  lassoOutput <- matrix(nrow=401, ncol=p)
  ListOfLambdas <- c()
  
  #Estimate lasso polynomial at every gridpoint
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    nonzeroWeights <- which(lassoWeights > 0)
    X <- buildFeature(currentPoint, p, x[nonzeroWeights])
    lassoFit <- cv.glmnet(X, y[nonzeroWeights], weights = lassoWeights[nonzeroWeights], maxit=10**7)
    ListOfLambdas <- c(ListOfLambdas, lassoFit$lambda.min)
    
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  
  #Smooth Lambdas output with loess, and refit with new lambda graph. This reduces variance
  smoothedLambas <- lowess(gridPoints, ListOfLambdas, f=1/10)$y
  smoothLassoOutput <- matrix(nrow=401, ncol=p)
  
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    nonzeroWeights <- which(lassoWeights > 0)
    X <- buildFeature(currentPoint, p, x[nonzeroWeights])
    lassoFit <- glmnet(X, y[nonzeroWeights], weights = lassoWeights[nonzeroWeights], standardize=TRUE, alpha=1, maxit=10**7)
    
    # Selects the desired coefficients with penalty lambda
    smoothLassoCoef <- coef(lassoFit, s=smoothedLambas[i])
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef) 
  }
  
  #Return list of smooth and unsmooth lasso, alongside the lambdas
  return(list(lasso=lassoOutput, smoothLasso=smoothLassoOutput,
              unsmoothedLambdas=ListOfLambdas, smoothedLambdas=smoothedLambas))
}




