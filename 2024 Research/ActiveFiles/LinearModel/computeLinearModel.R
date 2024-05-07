computeLinearModel <- function(observedX, observedY, bandwidth, evalPoints, degree){
  
  
  lmDerivatives <- matrix(nrow=length(evalPoints), ncol=degree+1)
  
  for(i in 1:length(evalPoints)){
    currentPoint <- evalPoints[i]
    
    #Weights and Feature matrix
    X <- buildFeatureDataFrame(currentPoint, degree, observedX)
    weights <- computeWeights(observedX, currentPoint, bandwidth)
    
    X$observedY <- observedY
    
    linearModel <- lm(observedY ~ .,data=X, weights=weights)
    print(linearModel$coefficients)
    
    lmDerivatives[i,] <- as.vector(linearModel$coefficients)
  }
  
  return(lmDerivatives)
  
}



# This function builds our feature matrix.
buildFeatureDataFrame <- function(midpoint, degree, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = degree)
  
  for(i in 1:degree){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(as.data.frame(X))
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