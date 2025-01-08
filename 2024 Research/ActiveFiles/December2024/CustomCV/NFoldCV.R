NfoldCV <-function(x, y, gridPoint, lambdas, weights, p, trueY, folds = 10){
  
  
  #Split data into n folds
  splitData <- split(1:length(x), sample(1:folds, length(x), replace=TRUE))
  #print(splitData)
  
  AVEWCV <- c()
  X <- buildFeature(gridPoint, p, x) #X does not change, we just have to leave out one row at a time
  
  
  for(lambda in lambdas){
    WMSE <- numeric(folds)
    
    for(i in 1:folds){
      
      currentLeftOut <- as.vector(splitData[[i]])
      #print(currentLeftOut)
      
      lassoFit <- glmnet::glmnet(X[-currentLeftOut,], y[-currentLeftOut], weights = weights[-currentLeftOut], maxit=10**7)
      lassoCoef <- as.vector(coef(lassoFit, s=lambda))
      evaluateModel <- buildPolynomial(x[currentLeftOut], lassoCoef, gridPoint)
      
      #print(trueY[-currentLeftOut])
      #print(class(trueY))
      #print(class(weights))

      #print(length(currentLeftOut))
      WMSE[i] <- sum(weights[currentLeftOut]*(trueY[currentLeftOut] - evaluateModel)**2)
    }
    #print(lambda)
    AVEWCV <- c(AVEWCV, mean(WMSE))
  }
  
  #Specify Best Lambda
  lambdaIndex <- which.min(AVEWCV)
  chosenLambda <- lambdas[lambdaIndex]
  #print(AVEWCV)
  #print(lambdaIndex)
  
  #Use lambda on full dataset
  lassoFit <- glmnet::glmnet(X, y, weights = weights, maxit=10**7)
  lassoCoef <- coef(lassoFit, s=lambda, exact=TRUE, x=X, y=y, weights=weights)
  
  
  return(list(lassoCoef, chosenLambda))
}