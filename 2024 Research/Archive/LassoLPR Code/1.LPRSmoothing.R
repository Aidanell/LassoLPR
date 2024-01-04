library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)

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


#Computes Nadaraya-Watson and Local Linear estimators of a specified derivative.
#Errors are also returned 
ComputeLocPolyMethods <- function(x, y, trueEquation, drv, bandwidth){
  
  NWFit <- locpoly(x, y, degree = drv, drv = drv, bandwidth = bandwidth)
  NWError <- findError(trueEquation, NWFit$x, NWFit$y)
  NWIntError <- findInteriorError(trueEquation, NWFit$x, NWFit$y)
  
  LLFit <- locpoly(x, y, degree = drv+1, drv = drv, bandwidth = bandwidth)
  LLError <- findError(trueEquation, LLFit$x, LLFit$y)
  LLIntError <- findInteriorError(trueEquation, LLFit$x, LLFit$y)
  
  ErrorData <- list(NWFit = NWFit, LLFit = LLFit,
                    NWError = NWError, NWIntError = NWIntError,
                    LLError = LLError, LLIntError = LLIntError)
  
  return(ErrorData)
}


#Performs either lasso or ridge regression depending on alpha value.
#Returns a vector of estimated derivatives
RegularizationLPR <- function(X, yValues, weights, alpha){
  RegFit <- cv.glmnet(X, yValues, weights = weights, standardize=TRUE, alpha=alpha)
  # Selects the desired coefficients with penalty lambda
  RegCoef <- coef(RegFit, s="lambda.min")
  return(as.vector(RegCoef))
}



#Main Method----
LPRSmoothing <- function(
        equation,   # The mathematical expression our data comes from
        left,       # Left endpoint of window
        right,      # Right endpoint of window
        numTerms,   # Number of Terms in the Taylor polynomial
        sigma,      # Noise for randomness in generated data
        numPoints,   # How many data points to be generated
        plot = TRUE, #If you want the data to be plotted
        seed = NULL, #Optional if you want to generate data with specific seed
        bandMulti = 8 #The factor the bandwidth returned from dpill will be multiplied by
        ){
  
  lassoData <- vector(mode='list', length=4) #Initializing Data structure to return at the end
  if(is.null(seed) == FALSE){set.seed(seed)} #Set seed if specified

  
  #Calculates all the derivatives needed for plotting
  derivList <- derivCalc(func = equation, numDeriv = 5)
  
  ##############Building Simulated Data##################
  
  
  x <- sort(runif(numPoints, min=left, max=right))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  exactY <- eval(equation)
  y <- exactY + noise
  
  #dpill returns estimated optimal bandwidth for NW and LL estimators.
  #Lasso requires a larger bandwidth and should be multiplied by a specified factor
  bandwidth <- dpill(x,y)*bandMulti 
  
  #Points at which we calculate our non-parametric estimators
  evalPoints <- seq(left, right, length.out=401)
  
  
  ##############Computing Smoothing Function#############
  
  #Matrix which will contain all estimated derivatives for each point
  LassoOutput <- RidgeOutput <- matrix(nrow=401, ncol=numTerms+1)
  
  for(i in 1:length(evalPoints)){
    currentPoint <- evalPoints[i]
    
    #Indicies of points within neighbourhood
    neighbourIndicies <- which(abs(x-currentPoint) < bandwidth)
    
    X <- buildFeature(currentPoint, numTerms, x[neighbourIndicies]) # Feature matrix
    weights <- dnorm(x[neighbourIndicies], mean = currentPoint, sd = bandwidth) #Weights
    
    LassoCoef <- RegularizationLPR(X, y[neighbourIndicies], weights, alpha=1) #Lasso-LPR
    RidgeCoef <- RegularizationLPR(X, y[neighbourIndicies], weights, alpha=0) #Ridge-LPR
    
    
    #Appends Lasso and Ridge outputs to their respective matrices
    LassoOutput[i, ] <- LassoCoef
    RidgeOutput[i, ] <- RidgeCoef
  }
  
  
  
  #################Computing Errors##################
  
  #ErrList will contain all errors of each method for the function
  #and the first 3 derivatives
  ErrList <- vector(mode='list', length=4)
  #for the function and its first 3 derivatives
  for(i in 1:4){
    
    #Computes LL and NW Errors
    ErrorData <- ComputeLocPolyMethods(x=x, y=y, trueEquation = derivList[i],
                                     drv = i-1, bandwidth = bandwidth/bandMulti)
    
    #Calculates and Appends Ridge/Lasso MSE
    ErrorData$LassoError <- findError(derivList[i], evalPoints, LassoOutput[,i])
    ErrorData$LassoIntError <- findInteriorError(derivList[i], evalPoints, LassoOutput[,i])
    ErrorData$RidgeError <- findError(derivList[i], evalPoints, RidgeOutput[,i])
    ErrorData$RidgeIntError <- findInteriorError(derivList[i], evalPoints, RidgeOutput[,i])
    
    ErrList[[i]] <- ErrorData
  }
  
  
  
  ##################Plotting#########################
  
  x <- evalPoints #for plotting exact functions 
  
  if(plot == TRUE){
    par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(2, 2, 2, 2)) #For spacing of legend
    
    plot(x, eval(derivList[1]), main="Regression", ylab='y', type='l', lwd=3)
    lines(evalPoints, LassoOutput[,1], col='red')
    lines(evalPoints, RidgeOutput[,1], col='orange')
    lines(ErrList[[1]]$NWFit, col='blue')
    lines(ErrList[[1]]$LLFit, col='green')
    
    plot(x, eval(derivList[2]), main="Regression Of 1st Derivative", ylab='y', type='l', lwd=3)
    lines(evalPoints, LassoOutput[,2], col='red')
    lines(evalPoints, RidgeOutput[,2], col='orange')
    lines(ErrList[[2]]$NWFit, col='blue')
    lines(ErrList[[2]]$LLFit, col='green')
    
    plot(x, eval(derivList[3]), main = "Regression Of 2nd Derivative", ylab='y', type='l', lwd=3)
    lines(evalPoints, LassoOutput[,3], col='red')
    lines(evalPoints, RidgeOutput[,3], col='orange')
    lines(ErrList[[3]]$NWFit, col='blue')
    lines(ErrList[[3]]$LLFit, col='green')
    
    plot(x, eval(derivList[4]), main = "Regression Of 3rd Derivative", ylab='y', type='l', lwd=3)
    lines(evalPoints, LassoOutput[,4], col='red')
    lines(evalPoints, RidgeOutput[,4], col='orange')
    lines(ErrList[[4]]$NWFit, col='blue')
    lines(ErrList[[4]]$LLFit, col='green')
    
    #Legend
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    plot_colors <- c("black", "red","blue", "green", "orange")
    legend(x = "bottom",inset = 0,
           legend = c("True Function", "Lasso-Smoothing", "Nadaraya-Watson", "Local Linear", "Ridge-Smoothing"), 
           col=plot_colors, lwd=5, cex=1, horiz = TRUE, xpd=TRUE)
  }
  
  
  ##############Returning Data################################
  #This is just formatting the data to be easy to access for the simulations
  for(i in 1:4){
    lassoFits <- list(LassoFit0 = list("x" = evalPoints, "y" = LassoOutput[,1]), 
                      LassoFit1 = list("x" = evalPoints, "y" = LassoOutput[,2]), 
                      LassoFit2 = list("x" = evalPoints, "y" = LassoOutput[,3]),
                      LassoFit3 = list("x" = evalPoints, "y" = LassoOutput[,4]))
    RidgeFits <- list(RidgeFit0 = list("x" = evalPoints, "y" = RidgeOutput[,1]), 
                      RidgeFit1 = list("x" = evalPoints, "y" = RidgeOutput[,2]), 
                      RidgeFit2 = list("x" = evalPoints, "y" = RidgeOutput[,3]),
                      RidgeFit3 = list("x" = evalPoints, "y" = RidgeOutput[,4]))
    
    
    data <- list("LassoFit" = lassoFits[[i]], "NWFit" = ErrList[[i]]$NWFit,
                 "LLFit" = ErrList[[i]]$LLFit, "RidgeFit" = RidgeFits[[i]],
                 "LassoError"=ErrList[[i]]$LassoError, "RidgeError"=ErrList[[i]]$RidgeError,
                 "NWError" = ErrList[[i]]$NWError, "LLError" = ErrList[[i]]$LLError,
                 "LassoIntError" = ErrList[[i]]$LassoIntError,
                 "RidgeIntError"=ErrList[[i]]$RidgeIntError,
                 "NWIntError"=ErrList[[i]]$NWIntError,
                 "LLIntError"=ErrList[[i]]$LLIntError)
    lassoData[[i]] <- data
  }
  return(lassoData)
}

#The three expressions that were tested
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))
sine <- expression(sin(5*pi*x))
bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))

#Test Run
testLasso <- LPRSmoothing(equation = peak,
                         left = 0,
                         right = 1,
                         numTerms = 10,
                         sigma = sqrt(0.5),
                         numPoints = 500)
testLasso





