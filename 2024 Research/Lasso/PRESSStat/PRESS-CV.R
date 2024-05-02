
##Running PRESS
#observedX: observe x
#observedY: observe y
#trueY: true y
#X: Feature matrix
#yhat: model output

PRESS.cv <- function(observedX, observedY, X, weights){
  
  #Fit lasso model at all lambdas. 
  fit = glmnet(X, observedY, weights=weights, standardize=TRUE, alpha=1, maxit=10**7)
  lambdas <- fit$lambda
  
  AllPressStats <- rep(0, length(lambdas))
  totalWeights <- sum(weights)
  
  #Hat matrix
  #print(det(t(X) %*% X))
  H <- X %*% MASS::ginv(t(X) %*% X) %*% t(X)
  Hdiag <- diag(H)

  #Calculate PRESS stat for each lambda
  for(i in 1:length(lambdas)){
    
    yhat <- predict(fit, newx=X, s=lambdas[i])
    
    PRESSVector <- ((yhat-observedY)/(1-Hdiag))**2 * (weights/totalWeights)
    PRESSStat <- sum(PRESSVector)
    AllPressStats[i] <- PRESSStat
  }
  #Getting best choice of lambda and returning corresponding coef
  minLambda <- lambdas[which.min(AllPressStats)]
  bestModelCoef <- as.vector(coef(fit, s=minLambda))

  PRESS.Data <- list(LassoCoef = bestModelCoef, minLambda = minLambda)
  return(PRESS.Data)
}

