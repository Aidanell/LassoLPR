




AidansCV <-function(x, y, gridPoint, lambdas, weights, p, trueY){
  
  
  #lambdas <- seq(2,0,length.out=20)
  AVEWCV <- numeric(20)
  X <- buildFeature(gridPoint, p, x) #X does not change, we just have to leave out one row at a time
  trueY <- unlist(trueY)
  
  
  
  for(lambda in lambdas){
    WMSE <- numeric(length(y))
    for(i in 1:length(y)){
      
      currentLeftOut <- y[i]
      
      lassoFit <- glmnet::glmnet(X[-i,], y[-i], weights = weights[-1], maxit=10**7)
      lassoCoef <- as.vector(coef(lassoFit, s=lambda))
      evaluateModel <- buildPolynomial(x[-1], lassoCoef, gridPoint)
      
      #print(length(trueY[-1]))
      #print(length(weights[-1]))
      #print(length(evaluateModel))
      
      WMSE[i] <- sum(weights[-1]*(trueY[-1] - evaluateModel)**2)
    }
    #print(lambda)
    AVEWCV <- mean(WMSE)
  }
  
  #Specify Best Lambda
  lambdaIndex <- which.min(AVEWCV)
  chosenLambda <- lambdas[lambdaIndex]
  print(chosenLambda)
  #Use lambda on full dataset
  lassoFit <- glmnet::glmnet(X, y, weights = weights, maxit=10**7)
  lassoCoef <- coef(lassoFit, s=lambda)
  
  
  return(list(lassoCoef, chosenLambda))
}



buildFeature <- function(midpoint, p, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = p)
  for(i in 1:p){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}


buildPolynomial <- function(x, coefs, gridPoint){
  output <- coefs[1]
  for(i in 2:length(coefs)){
    output <- output + coefs[i]*(x-gridPoint)^i / factorial(i)
  }
  return(output)
}


#Function of the Epanechnikov kernel
epan <- function(x, bandwidth){
  ifelse(abs(x) <= 1, 3/4 * (1 - x^2), 0)
}


#Computes the scaled kernel weights for every input x
#Works for Epanechnikov or normal kernels
computeWeights <- function(x, x0, bandwidth, kernel='epak'){
  weights <- numeric(length(x))
  
  for(i in 1:length(x)){
    diff <- abs(x[i]-x0) / bandwidth
    if(kernel == 'epak'){
      weights[i] <- epan(diff)/bandwidth
    }else{
      weights[i] <- dnorm(diff)/bandwidth
    }
  }
  return(weights)
}

curve(buildPolynomial(x, coefs=c(1,2,3), gridPoint=2), from=0, to=4)

