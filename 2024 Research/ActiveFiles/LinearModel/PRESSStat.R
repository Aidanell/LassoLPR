computePRESSLasso <- function(observedX, observedY, bandwidth, evalPoints, degree){
  
  
  PRESSDerivatives <- matrix(nrow=length(evalPoints), ncol=degree+1)
  lambdas <- c()
  
  for(i in 1:length(evalPoints)){
    currentPoint <- evalPoints[i]
    
    #Weights and Feature matrix
    X <- buildFeature(currentPoint, degree, observedX)
    weights <- computeWeights(observedX, currentPoint, bandwidth)
    
    PRESSCV <- PRESS.cv(observedX, observedY, X, weights)
    
    lambdas <- c(lambdas, PRESSCV$minLambda)
    PRESSDerivatives[i,] <- PRESSCV$LassoCoef
  }
  
  PRESSLASSOData <- list(PRESSDerivatives = PRESSDerivatives, lambdas = lambdas)
  return(PRESSLASSOData)
  
}


# This function builds our feature matrix.
buildFeature <- function(midpoint, degree, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = degree)
  
  for(i in 1:degree){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}



#Computes the scaled kernel weights for every input x
computeWeights <- function(x, midpoint, bandwidth){
  output <- vector(length=length(x))
  for(i in 1:length(x)){
    
    diff <- abs(x[i]-midpoint) / bandwidth
    if(diff > 1){output[i] <- 0}
    else{output[i] <- epan(diff)/ (bandwidth)}
    
  }
  return(output)
}


#Function of the Epanechnikov kernel
epan <- function(x, bandwidth){
  return(3/4 * (1-x**2))
}


#Calculates true derivatives needed to compare estimates of derivatives with real
derivCalc <- function(func, numDeriv){
  drvList <- c(func)
  nextDrv <- func
  for(i in 1:numDeriv){
    nextDrv <- D(nextDrv, 'x')
    drvList <- c(drvList, nextDrv)
  }
  return(drvList)
}
