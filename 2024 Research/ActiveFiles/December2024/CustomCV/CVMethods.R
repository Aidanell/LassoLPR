




AidansCV <-function(x, y, gridPoint, lambdas, weights, p, trueFunctionGridPoints){
  
  
  lambdas <- seq(2,0,length.out=20)
  AVE-WCV <- numeric(20)
  X <- buildFeature(gridPoint, p, x) #X does not change, we just have to leave out one row at a time
  
  for(lambda in lambdas){
    
    for(i in 1:length(y)){
      WMSE <- numeric(length(y))
      
      currentLeftOut <- y[i]
      
      lassoFit <- glmnet::glmnet(X[-i,], y[-i], weights = lassoWeights[-1], maxit=10**7)
      lassoCoef <- coef(lassoFit, s=lambda, exact=TRUE)
      
      WMSE[i] <- sum(lassoWeights[-1]*(trueFunctionGridPoints - buildPolynomial(x[-1], lassoCoef, gridPoint))**2)
    }
    AVE-WCV <- mean(CV-WMSE)
  }
  
  lambdaIndex <- which.min(AVE-WCV)
  return(lambdas[lambdaIndex])
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


curve(buildPolynomial(x, coefs=c(1,2,3), gridPoint=2), from=0, to=4)

