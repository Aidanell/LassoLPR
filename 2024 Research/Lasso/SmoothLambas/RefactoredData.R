library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)



#Main Method----
RunLPRSimulation<- function(
    equation,   # The mathematical expression our data comes from
    left,       # Left endpoint of window
    right,      # Right endpoint of window
    numTerms,   # Number of Terms in the Taylor polynomial
    sigma,      # Noise for randomness in generated data
    numPoints,   # How many data points to be generated
    plot = TRUE, #If you want the data to be plotted
    seed = NULL, #Optional if you want to generate data with specific seed
    lassoBandMulti = 8, #The factor the bandwidth returned from dpill will be multiplied by for lasso
    ridgeBandMulti = 5 #The factor the bandwidth returned from dpill will be multiplied by for ridge,
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
  
  computeSmoothLasso <- function(x,y ,lassoBandwidth, ridgeBandwidth, 
                                           evalPoints, numTerms, lambdas){
    
    lassoOutput <- matrix(nrow=401, ncol=numTerms+1)
    ridgeOutput <- matrix(nrow=401, ncol=numTerms+1)
    ListOfLambdas <- c()
    
    
    for(i in 1:length(evalPoints)){
      currentPoint <- evalPoints[i]
      
      X <- buildFeature(currentPoint, numTerms, x)
      
      lassoWeights <- computeWeights(x, currentPoint, lassoBandwidth)
      ridgeWeights <- computeWeights(x, currentPoint, ridgeBandwidth)
      
      lassoFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
      ridgeFit <- glmnet(X, y, weights = ridgeWeights, standardize=TRUE, alpha=0)
      
      ListOfLambdas <- c(ListOfLambdas, lambdas[i])
      
      # Selects the desired coefficients with penalty lambda
      lassoCoef <- coef(lassoFit, s=lambdas[i])
      ridgeCoef <- coef(ridgeFit, s=lambdas[i])
      
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
  derivList <- derivCalc(func = equation, numDeriv = numTerms)
  
  ##############Building Simulated Data##################
  
  x <- sort(runif(numPoints, min=left, max=right))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  exactY <- eval(equation)
  y <- exactY + noise
  
  #dpill returns estimated optimal bandwidth for locpoly estimators.
  #Lasso requires a larger bandwidth and should be multiplied by a specified factor
  locpolyBandwidth <- dpill(x,y)
  lassoBandwidth <- locpolyBandwidth*lassoBandMulti
  ridgeBandwidth <- locpolyBandwidth*ridgeBandMulti
  
  #Points at which we calculate our non-parametric estimators
  evalPoints <- seq(left, right, length.out=401)
  
  ##############Computing Smoothing Function#############
  
  regLPRMethods <- computeRegularizationMethods(x, y, lassoBandwidth, ridgeBandwidth,
                                                evalPoints, numTerms)
  lassoOutput <- regLPRMethods[[1]]
  ridgeOutput <- regLPRMethods[[2]]
  listOfLambdas <- regLPRMethods[[3]]
  
  ##########Computing Smoothed Lambdas Lasso############
  smoothedLambas <- lowess(evalPoints, listOfLambdas, f=1/10)$y
  smoothLambaMethod <- computeSmoothLasso(x, y, lassoBandwidth, ridgeBandwidth,
                                          evalPoints, numTerms, smoothedLambas)
  
  smoothLambdaOutput <- smoothLambaMethod[[1]]
  smoothedLambdas <- smoothLambaMethod[[3]]
  
  #################Computing Errors##################
  
  #ErrList will contain all errors of each method for the function and the first 3 derivatives
  ErrorTable <- matrix(nrow=0, ncol=8)
  colnames(ErrorTable) <- c("lassoError", "lassoIntError", "ridgeError", "ridgeIntError",
                            "locpolyError", "locpolyIntError", "smoothLamError", "smoothLamIntError")
  locpolyOutput <- matrix(nrow=401, ncol=numTerms+1)
  
  for(i in 1:11){
    
    #Computes locpoly estimation and it's errors
    locpolyFit <- KernSmooth::locpoly(x, y, degree = i-1, drv = i-1, bandwidth = locpolyBandwidth , gridsize = 401)
    
    #Calculates and Appends Ridge/Lasso MSE to datatable
    lassoError <- findError(derivList[i], evalPoints, lassoOutput[,i])
    lassoIntError <- findInteriorError(derivList[i], evalPoints, lassoOutput[,i])
    ridgeError <- findError(derivList[i], evalPoints, ridgeOutput[,i])
    ridgeIntError <- findInteriorError(derivList[i], evalPoints, ridgeOutput[,i])
    locpolyError <- findError(derivList[i], locpolyFit$x, locpolyFit$y)
    locpolyIntError <- findInteriorError(derivList[i], locpolyFit$x, locpolyFit$y)
    smoothLambdaError <- findError(derivList[i], evalPoints, smoothLambdaOutput[,i])
    smoothLambdaIntError <- findInteriorError(derivList[i], evalPoints, smoothLambdaOutput[,i])
    
    ErrorTable <- rbind(ErrorTable, c(lassoError, lassoIntError, ridgeError,
                                      ridgeIntError, locpolyError, locpolyIntError,
                                      smoothLambdaError, smoothLambdaIntError))
    locpolyOutput[,i] <- locpolyFit$y
    
  }
  xpoints <- x
  x <- evalPoints #for plotting exact functions 
  
  ##############Returning Data################################
    
  SimulatedData <- list("x" = xpoints, "y" = y)
    
  lassoObject <- list("simulatedData"=SimulatedData, "lassoOutput"=lassoOutput,
                      "ridgeOutput"=ridgeOutput, "locpolyOutput"=locpolyOutput,
                      "smoothLassoOutput"=smoothLambdaOutput,"errorTable"=ErrorTable,
                      "evalPoints" = x,
                      "lambdas" = list("unsmoothed"=listOfLambdas, "smoothed"=smoothedLambdas),
                      "derivFunctions"=derivList, "seed"=seed)
  
  return(lassoObject)
}


plotData <- function(simObject, deriv=0){
  
  x <- simObject$evalPoints
  trueFunction <- simObject$derivFunctions[deriv+1]
  y <- eval(trueFunction)
  
  if(length(y) == 1){ #For when trueFunction is constant
    y <- rep(y, 401)
  }
  
  lassoY <- simObject$lassoOutput[,deriv+1]
  ridgeY <- simObject$ridgeOutput[,deriv+1]
  locpolyY <- simObject$locpolyOutput[,deriv+1]
  smoothLassoY <- simObject$smoothLassoOutput[,deriv+1]
  
  combinedOutputs <- c(lassoY, ridgeY, smoothLassoY)
  plotLimits <- c(min(combinedOutputs), max(combinedOutputs))
  plotTitle <- paste("Derivative =", deriv)
  plot(x,y, type='l', lwd=3, ylim=plotLimits, main=plotTitle)
  lines(x, lassoY, col='red', lwd=2)
  lines(x, ridgeY, col='green', lwd=2)
  lines(x, locpolyY, col='blue', lwd=2)
  lines(x, smoothLassoY, col='orange', lwd=2)
  
}

#Running the script----

#The three expressions that were tested
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))
sine <- expression(sin(5*pi*x))
bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))



# poly <- expression(x**3 - 1.5*x**2 + 0.7)
# left = -0.6
# right = 1.4

#Run a Lasso vs Ridge vs Locpoly simulation
testLassoEpan <- RunLPRSimulation(equation = peak,
                                  left = 0,
                                  right = 1,
                                  numTerms = 10,
                                  sigma = sqrt(0.5),
                                  numPoints = 500,
                                  seed=50)

testLassoEpan
par(mfrow=c(2,2))
for(i in 0:3){
  plotData(testLassoEpan, deriv=i)
}

