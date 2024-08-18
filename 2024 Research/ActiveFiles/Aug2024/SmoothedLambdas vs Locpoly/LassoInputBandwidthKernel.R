


LassoSmoothed <- function(x,y, bandwidth, epanKernal=TRUE, p){
  
  evalPoints <- seq(0,1, length.out = 401)
  lassoOutput <- matrix(nrow=401, ncol=p+1)
  unsmoothedLambdas <- c()
  
  for(i in 1:length(evalPoints)){
    currentPoint <- evalPoints[i]
    
    X <- buildFeature(currentPoint, numTerms, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth, epan=epanKernal)
    lassoFit <- cv.glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
    unsmoothedLambdas <-  c(unsmoothedLambdas, lassoFit$lambda.min)
    
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  
  
  smoothedLambas <- lowess(evalPoints, unsmoothedLambdas, f=1/10)$y
  smoothLassoOutput <- matrix(nrow=401, ncol=numTerms+1)
  
  for(i in 1:length(evalPoints)){
    currentPoint <- evalPoints[i]
    
    X <- buildFeature(currentPoint, numTerms, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth, epan=epanKernal)
    lassoFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
    
    # Selects the desired coefficients with penalty lambda
    smoothLassoCoef <- coef(lassoFit, s=smoothedLambas[i])
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef) 
  }
  
  return(list(lassoOutput, smoothLassoOutput, unsmoothedLambdas, smoothedLambas))
}






