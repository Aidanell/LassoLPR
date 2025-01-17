NfoldCV <-function(x, y, gridPoint, lambdas, weights, p, trueY, folds = 10){
  
  
  #Split data into n folds
  splitData <- split(1:length(x), sample(1:folds, length(x), replace=TRUE))
  
  AVEWCV <- c()
  X <- buildFeature(gridPoint, p, x) #X does not change, we just have to leave out one row at a time
  
  
  for(lambda in lambdas){
    WMSE <- numeric(folds)
    
    for(i in 1:folds){
      
      currentLeftOut <- as.vector(splitData[[i]])
      #print(currentLeftOut)
      
      lassoFit <- glmnet::glmnet(X[-currentLeftOut,], y[-currentLeftOut], weights = weights[-currentLeftOut], maxit=10**7)
      lassoCoef <- as.vector(coef(lassoFit, s=lambda))
      #Evaluate model with points that were not trained with
      evaluateModel <- buildPolynomial(x[currentLeftOut], lassoCoef, gridPoint)
      
      #Get weighted MSE (with more weight towards points closer to gridPoint)
      WMSE[i] <- sum(weights[currentLeftOut]*(trueY[currentLeftOut] - evaluateModel)**2)
    }
    #Add CVMSE to list
    AVEWCV <- c(AVEWCV, mean(WMSE))
  }
  
  #Specify Best Lambda
  lambdaIndex <- which.min(AVEWCV)
  chosenLambda <- lambdas[lambdaIndex]
  
  #Use lambda on full dataset
  lassoFit <- glmnet::glmnet(X, y, weights = weights, maxit=10**7)
  lassoCoef <- coef(lassoFit, s=lambda, exact=TRUE, x=X, y=y, weights=weights)
  
  
  return(list(lassoCoef, chosenLambda))
}