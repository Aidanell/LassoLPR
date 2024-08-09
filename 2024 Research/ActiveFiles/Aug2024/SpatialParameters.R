##This file is intended to calculate an optimal data-adaptive ridge parameter
#with the algorithim proposed in Ch 2.1 of Siefert and Gaser (1998)

#We need sigma, s2, s1, s0, s2*, delta, and first derivative at xTilde


calcSj <- function(j, x, xTilde, kernelWeights){
  
  Sj <- sum(kernelWeights * (x - xTilde)**j)
  return(Sj)
}

calcSjStar <- function(j, x, xTilde, kernelWeights){
  
  SjStar <- sum(kernelWeights**2 * (x - xTilde)**j)
  return(SjStar)
}


peakDerivative <- function(x){
  return(-4000*(x-0.5)*exp(-4000*(x-0.5)**2) - 5)
}


GetSpatialRidgeList <- function(x, evalPoints, bandwidth, XTildeTracking, firstDeriv){
  
  SpatialRidgeList <- c()
  for(i in 1:length(evalPoints)){
    x0 <- evalPoints[i]
    lassoWeights <- computeWeights(x, x0, lassoBandwidth)
    xTilde <- XTildeTracking[i]
    rXtilde <- firstDeriv[i]
    
    nextRbar <- singleSpatialRidge(x, xTilde, lassoWeights, rXtilde, x0)
    
    SpatialRidgeList <- c(SpatialRidgeList, nextRbar)
    
  }
  
  return(SpatialRidgeList)
}

singleSpatialRidge <- function(x, xTilde, weights, rXtilde, x0){
  
  sigma <- sqrt(0.5) ####Change to estimate with data when I figure out the right method
  delta <- x0 - xTilde
  
  s0 <- calcSj(0, x, xTilde, weights)
  s1 <- calcSj(1, x, xTilde, weights)
  s2 <- calcSj(2, x, xTilde, weights)
  s1Star <- calcSjStar(1, x, xTilde, weights)
  s2Star <- calcSjStar(2, x, xTilde, weights)

  
  
  #Now use formula to get RbarOpt
  
  numerator <- (rXtilde * delta)**2 - sigma**2 * (delta * s1Star)/(s0 * s2)
  denomenator <- (rXtilde * delta)**2 + sigma**2 * (delta**2 * s2Star)/s2**2
  
  RbarOpt <- numerator / denomenator
  
  Rbar <- s2 * (1-RbarOpt)/RbarOpt
  
  if(RbarOpt > 0 && RbarOpt < 1){
    return(RbarOpt)
  }else{
    return(Rbar) 
  }
}
