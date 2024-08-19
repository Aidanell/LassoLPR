


buildLassoMatricies <- function(x,y,z,n,p, Xbandwidth, Ybandwidth, evalPointsX, evalPointsY, gridLength){
  
  ListOfLassoMatricies <- vector(mode='list', length=ncolFunction(p))
  unSmoothedLambdas <- rep(0,length(evalPointsX))
  
  for(i in 1:ncolFunction(p)){ListOfLassoMatricies[[i]] <- rep(0,length(evalPointsX))}
  
  for(i in 1:length(evalPointsX)){
    #current coordinates of centered Taylor polynomial
    currentX <- evalPointsX[i]
    currentY <- evalPointsY[i]
    
    #Build feature matrix and weights
    feature <- build2DFeature(currentX, currentY, x, y, p)
    #print(feature)
    xWeights <- computeWeights(x, currentX, Xbandwidth)
    yWeights <- computeWeights(y, currentY, Ybandwidth)
    weights <- xWeights * yWeights
    
    #Fit lasso and find best lambda
    lassoFit <- cv.glmnet(feature, z, weights = weights, standardize=TRUE, maxit=1000000)
    lassoCoef <- as.vector(coef(lassoFit, s="lambda.min"))
    
    #Organize all the derivative estimates
    for(k in 1:ncolFunction(p)){
      ListOfLassoMatricies[[k]][i] <- lassoCoef[k]
    }
    unSmoothedLambdas[i] <- lassoFit$lambda.min
    
    if(i %% gridLength == 0){print(paste(i/gridLength,"/",gridLength))}
  }
  
  return(list(ListOfLassoMatricies, unSmoothedLambdas))
}


build2DFeature <- function(centerX, centerY, xValues, yValues, p){
  #This algorithim incrementally develops the coefficients needed
  #to build taylor polynomials up to a desired p. Constant column not included
  #as it's built with glmnet
  featureMatrix <- matrix(nrow=length(xValues), ncol=ncolFunction(p)-1)
  
  currentP <- 1
  currentCol <- 1
  
  while(currentP <= p){
    numRows <- currentP
    
    for(t in numRows:0){
      numer <- (xValues - centerX)**2 * (yValues - centerY)**(numRows-t)
      denom <- factorial(t)*factorial(numRows-t)
      
      featureMatrix[,currentCol] <- numer/denom
      currentCol  <- currentCol + 1
    }
    currentP <- currentP + 1
  }

  return(featureMatrix)
}


computeWeights <- function(x, midpoint, bandwidth){
  
  epan <- function(x, bandwidth){return(3/4 * (1-x**2))}
  
  output <- vector(length=length(x))
  for(i in 1:length(x)){
    
    diff <- abs(x[i]-midpoint) / bandwidth
    if(diff > 1){output[i] <- 0}
    else{output[i] <- epan(diff)/ (bandwidth)}
    
  }
  return(output)
}


smoothLassoComputation<- function(x,y,z,n,p, smoothedLambdas, Xbandwidth, Ybandwidth, evalPointsX, evalPointsY, gridLength){
  
  smoothedLassoMatricies <- vector(mode='list', length=ncolFunction(p))
  
  for(i in 1:ncolFunction(p)){smoothedLassoMatricies[[i]] <- rep(0,length(evalPointsX))}
  
  for(i in 1:length(evalPointsX)){
    #current coordinates of centered Taylor polynomial
    currentX <- evalPointsX[i]
    currentY <- evalPointsY[i]
    
    #Build feature matrix and weights
    feature <- build2DFeature(currentX, currentY, x, y, p)
    #print(feature)
    xWeights <- computeWeights(x, currentX, Xbandwidth)
    yWeights <- computeWeights(y, currentY, Ybandwidth)
    weights <- xWeights * yWeights
    
    #Fit lasso and find best lambda
    lassoFit <- glmnet(feature, z, weights = weights, standardize=TRUE, maxit=1000000)
    lassoCoef <- as.vector(coef(lassoFit, s=smoothedLambdas[i]))
    
    #Organize all the derivative estimates
    for(k in 1:ncolFunction(p)){
      smoothedLassoMatricies[[k]][i] <- lassoCoef[k]
    }
    
    if(i %% gridLength == 0){print(paste(i/gridLength,"/",gridLength))}
  }
  
  return(smoothedLassoMatricies)
}


getTrueFunctionData <- function(xylim, gridLength){
  #Functions for comparing derivatives----
  #NOTE: ONLY FOR PEAK FUNCTION'S FIRST THREE DERIVATIVES
  #FUNCTIONS MUST BE MANUALLY CHANGED TO COMPARE OTHER FUNCTIONS
  f <- function(x,y){return(2-5*x-5*y +5*exp(-20*(x-0.5)**2 -20*(y-0.5)**2))}
  
  fx <- function(x,y){return(-200*(x-0.5)*exp(-20*(x-0.5)**2 -20*(y-0.5)**2) - 5)}
  fxx <- function(x,y){return((8000*(x**2) - 8000*x + 1800)*exp(-20*(x-0.5)**2 - 20*(y-0.5)**2))}
  fxxx <- function(x,y){return(-1*(320000*(x**3) - 480000*(x**2) + 216000*x - 28000)*exp(-20*(x-0.5)**2 -20*(y-0.5)**2))}
  
  fy <- function(x,y){return(-200*(y-0.5)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2) - 5)}
  fyy <- function(x,y){return((8000*(y**2) - 8000*y + 1800)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2))}
  fyyy <- function(x,y){return(-1*(320000*(y**3) - 480000*(y**2) + 216000*y - 28000)*exp(-20*(y-0.5)**2 -20*(x-0.5)**2))}
  
  fxy <- function(x,y){return(2000*(2*x - 1)*(2*y - 1)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2))}
  fxxy <- function(x,y){return(-4000*(40*(x**2) - 40*x +9)*(2*y - 1)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2))}
  fxyy <- function(x,y){return(-4000*(40*(y**2) - 40*y +9)*(2*x - 1)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2))}
  
  #Defining matricies of true data for the function and all its derivatives----
  trueX <- trueY <- seq(xylim[1], xylim[2], length.out = gridLength)
  trueResult <- vector(mode='list', length=10)
  trueResult[[1]] <- outer(trueX, trueY, f)
  trueResult[[2]] <- outer(trueX, trueY, fx)
  trueResult[[3]] <- outer(trueX, trueY, fy)
  trueResult[[4]] <- outer(trueX, trueY, fxx)
  trueResult[[5]] <- outer(trueX, trueY, fxy)
  trueResult[[6]] <- outer(trueX, trueY, fyy)
  trueResult[[7]] <- outer(trueX, trueY, fxxx)
  trueResult[[8]] <- outer(trueX, trueY, fxxy)
  trueResult[[9]] <- outer(trueX, trueY, fxyy)
  trueResult[[10]] <- outer(trueX, trueY, fyyy)
  
  return(trueResult)
  
}

#How many columns needed for lasso Taylor polynomials of degree p
ncolFunction <- function(p){return(1+p+sum(1:p))}






