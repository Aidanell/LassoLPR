
#This file contains helper methods used in cv.glmnet calculations

# This function builds our feature matrix.
buildFeature <- function(midpoint, numTerms, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = numTerms)
  
  for(i in 1:numTerms){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
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
