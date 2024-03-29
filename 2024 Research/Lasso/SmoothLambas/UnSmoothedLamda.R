UnsmoothedLambdasLPR <- function(
    equation,   # The mathematical expression our data comes from
    left,       # Left endpoint of window
    right,      # Right endpoint of window
    numTerms,   # Number of Terms in the Taylor polynomial
    sigma,      # Noise for randomness in generated data
    numPoints,   # How many data points to be generated
    plot = TRUE, #If you want the data to be plotted
    seed = NULL, #Optional if you want to generate data with specific seed
    lassoBandMulti = 8, #The factor the bandwidth returned from dpill will be multiplied by for lasso
    ridgeBandMulti = 5 #The factor the bandwidth returned from dpill will be multiplied by for ridge
){
  
  #Helper Methods----
  
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
  
  
  #Computes MSE
  findError <- function(equation, estimateX, estimateY){
    x <- estimateX
    trueY <- eval(equation)
    MSE <- mean((estimateY - trueY)^2)
    return(MSE)
  }
  
  
  
  #Excludes the edges of the function in the MSE calculation
  findInteriorError <- function(equation, estimateX, estimateY, percentageCut=0.05){
    leftCutoff = as.integer(length(estimateX) * percentageCut)
    rightCutoff = as.integer(length(estimateY) - leftCutoff)
    
    x <- estimateX[leftCutoff : rightCutoff]
    y <- estimateY[leftCutoff : rightCutoff]
    trueY <- eval(equation)
    MSE <- mean((y - trueY)^2)
    return(MSE)
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
  
  
  #Computes LocPoly estimators of a specified derivative.
  #Errors are also returned 
  ComputeLocPolyMethods <- function(x, y, trueEquation, drv, bandwidth){
    
    locpolyFit <- KernSmooth::locpoly(x, y, degree = drv, drv = drv, bandwidth = bandwidth, gridsize = 401)
    locpolyError <- findError(trueEquation, locpolyFit$x, locpolyFit$y)
    locpolyIntError <- findInteriorError(trueEquation, locpolyFit$x, locpolyFit$y)
    
    ErrorData <- list(locpolyFit = locpolyFit, locpolyError = locpolyError,
                      locpolyIntError = locpolyIntError)
    return(ErrorData)
  }
  
  
  #Function of the Epanechnikov kernel
  epan <- function(x, bandwidth){
    return(3/4 * (1-x**2))
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
  
  
  computeRegularizationMethods <- function(x,y ,lassoBandwidth, ridgeBandwidth, 
                                           evalPoints, numTerms){
    
    lassoOutput <- matrix(nrow=401, ncol=numTerms+1)
    ridgeOutput <- matrix(nrow=401, ncol=numTerms+1)
    ListOfLambdas <- c()
    
    
    for(i in 1:length(evalPoints)){
      currentPoint <- evalPoints[i]
      
      X <- buildFeature(currentPoint, numTerms, x)
      
      lassoWeights <- computeWeights(x, currentPoint, lassoBandwidth)
      ridgeWeights <- computeWeights(x, currentPoint, ridgeBandwidth)
      
      lassoFit <- cv.glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
      ridgeFit <- cv.glmnet(X, y, weights = ridgeWeights, standardize=TRUE, alpha=0)
      
      ListOfLambdas <- c(ListOfLambdas, lassoFit$lambda.min)
      
      # Selects the desired coefficients with penalty lambda
      lassoCoef <- coef(lassoFit, s="lambda.min")
      ridgeCoef <- coef(ridgeFit, s="lambda.min")
      
      #The i'th row contains all derivative estimation at that specificed evalpoint
      lassoOutput[i, ] <- as.vector(lassoCoef) 
      ridgeOutput[i, ] <- as.vector(ridgeCoef)
    }
    return(list(lassoOutput, ridgeOutput, ListOfLambdas))
  }
  
  ####Setting Seed#####
  lassoData <- vector(mode='list', length=6) #Initializing Data structure to return at the end
  if(is.null(seed) == FALSE){set.seed(seed)} #Set seed if specified
  else{
      seed <- runif(1, min=0, max=9000000)
      set.seed(seed)
    }
  
  
  #Calculates all the derivatives needed for plotting exact functions and calculating errors
  derivList <- derivCalc(func = equation, numDeriv = 3)
  
  ##############Building Simulated Data##################
  
  x <- sort(runif(numPoints, min=left, max=right))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  exactY <- eval(equation)
  y <- exactY + noise
  
  #dpill returns estimated optimal bandwidth for locpoly estimators.
  #Lasso requires a larger bandwidth and should be multiplied by a specified factor
  locpolyBandwdith <- dpill(x,y)
  lassoBandwidth <- locpolyBandwdith*lassoBandMulti
  ridgeBandwidth <- locpolyBandwdith*ridgeBandMulti
  
  #Points at which we calculate our non-parametric estimators
  evalPoints <- seq(left, right, length.out=401)
  
  ##############Computing Smoothing Function#############
  
  regLPRMethods <- computeRegularizationMethods(x, y, lassoBandwidth, ridgeBandwidth,
                                                evalPoints, numTerms)
  lassoOutput <- regLPRMethods[[1]]
  ridgeOutput <- regLPRMethods[[2]]
  listOfLambdas <- regLPRMethods[[3]]
  
  #################Computing Errors##################
  
  #ErrList will contain all errors of each method for the function and the first 3 derivatives
  ErrList <- vector(mode='list', length=4)
  
  for(i in 1:4){
    
    #Computes locpoly estimation and it's errors
    ErrorData <- ComputeLocPolyMethods(x=x, y=y, trueEquation = derivList[i],
                                       drv = i-1, bandwidth = locpolyBandwdith)
    
    #Calculates and Appends Ridge/Lasso MSE to datatable
    ErrorData$LassoError <- findError(derivList[i], evalPoints, lassoOutput[,i])
    ErrorData$LassoIntError <- findInteriorError(derivList[i], evalPoints, lassoOutput[,i])
    ErrorData$RidgeError <- findError(derivList[i], evalPoints, ridgeOutput[,i])
    ErrorData$RidgeIntError <- findInteriorError(derivList[i], evalPoints, ridgeOutput[,i])
    
    ErrList[[i]] <- ErrorData
  }
  
  ##############Returning Data################################
  xpoints <- x
  x <- evalPoints #for plotting exact functions 
  
  #This is just formatting the data to be easy to access for the simulations
  lassoFits <- list(LassoFit0 = list("x" = evalPoints, "y" = lassoOutput[,1]), LassoFit1 = list("x" = evalPoints, "y" = lassoOutput[,2]), 
                    LassoFit2 = list("x" = evalPoints, "y" = lassoOutput[,3]), LassoFit3 = list("x" = evalPoints, "y" = lassoOutput[,4]))
  RidgeFits <- list(RidgeFit0 = list("x" = evalPoints, "y" = ridgeOutput[,1]), RidgeFit1 = list("x" = evalPoints, "y" = ridgeOutput[,2]), 
                    RidgeFit2 = list("x" = evalPoints, "y" = ridgeOutput[,3]), RidgeFit3 = list("x" = evalPoints, "y" = ridgeOutput[,4]))
  
  for(i in 1:4){
    data <- list("LassoFit" = lassoFits[[i]], "LocpolyFit" = ErrList[[i]]$locpolyFit, "RidgeFit" = RidgeFits[[i]],
                 "LassoError"=ErrList[[i]]$LassoError, "RidgeError"=ErrList[[i]]$RidgeError,
                 "LocpolyError" = ErrList[[i]]$locpolyError, "LassoIntError" = ErrList[[i]]$LassoIntError, 
                 "RidgeIntError"=ErrList[[i]]$RidgeIntError, "LocpolyIntError"=ErrList[[i]]$locpolyIntError,
                 "x" = xpoints, "y" = y, "DerivFunc" = derivList[i], "EvalPoints" = x)
    lassoData[[i]] <- data
  }
  
  lassoData[[5]] <- listOfLambdas #For plotting lambas
  lassoData[[6]] <- seed
  return(lassoData)
}
