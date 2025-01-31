NFoldCV_V2 <-function(x, y, gridPoint, lambdas, weights, p, trueY, folds = 10){
  
  
  #Split data into n folds
  splitData <- split(1:length(x), sample(1:folds, length(x), replace=TRUE))
  
  #Demonstrates how weights are dispoportionate
  #totalWeight <- numeric(folds)
  #for(i in 1:folds){
  #  totalWeight[i] <- sum(weights[splitData[[i]]])
  #}
  #barplot(totalWeight)
  
  AVEWCV <- c()
  X <- buildFeature(gridPoint, p, x) #X does not change, we just have to leave out one row at a time
  
  
  for(lambda in lambdas){
    WMSE <- numeric(folds)
    
    for(i in 1:folds){
      
      currentLeftOut <- as.vector(splitData[[i]])
      
      #Fit model without one of the folds
      lassoFit <- glmnet::glmnet(X[-currentLeftOut,], y[-currentLeftOut], weights = weights[-currentLeftOut], maxit=10**7)
      lassoCoef <- as.vector(coef(lassoFit, s=lambda))
      #plot(lassoFit)
      
      
      #We only want to evaluate with points that are close to the gridPoint (and in the fold)
      #windowSize <- 0.1
      #pointstoEval <- which(abs(x[currentLeftOut] - gridPoint) < 0.1 )
      
      
      #Evaluate model with points that were not trained with
      #evaluateModel <- buildPolynomial(x[currentLeftOut], lassoCoef, gridPoint)
      evaluateModel <- predict(lassoFit, newx=X[currentLeftOut,], s=lambda)
      
      scaledWeights <- weights[currentLeftOut] / sum(weights[currentLeftOut])
      #Get weighted MSE (with more weight towards points closer to gridPoint)
      WMSE[i] <- sum(scaledWeights*(trueY[currentLeftOut] - evaluateModel)**2)
    }
    #Add CVMSE to list
    AVEWCV <- c(AVEWCV, mean(WMSE))
  }
  
  #Specify Best Lambda
  lambdaIndex <- which.min(AVEWCV)
  chosenLambda <- lambdas[lambdaIndex]
  
  #plot(1:length(lambdas),lambdas)
  #points(lambdaIndex, chosenLambda, col='red', cex=3, pch=16)
  
  #Use lambda on full dataset
  lassoFit <- glmnet::glmnet(X, y, weights = weights, maxit=10**7)
  #lassoCoef <- coef(lassoFit, s=lambda, exact=TRUE, x=X, y=y, weights=weights)
  lassoCoef <- coef(lassoFit, s=lambda)
  
  #plot(lassoFit, xvar='dev', label=TRUE)
  
  
  return(list(lassoCoef, chosenLambda))
}