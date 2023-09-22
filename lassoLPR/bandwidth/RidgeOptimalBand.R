library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)
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

computeRidgeOutput <- function(x, y, bandwidth, evalPoints, numTerms){
  
  ridgeOutput <- matrix(nrow=401, ncol=numTerms+1)
  for(i in 1:length(evalPoints)){
    currentPoint <- evalPoints[i]
    
    X <- buildFeature(currentPoint, numTerms, x)
    weights <- computeWeights(x, currentPoint, bandwidth)
    
    ridgeFit <- cv.glmnet(X, y, weights = weights, standardize=TRUE, alpha=0)
    # Selects the desired coefficients with penalty lambda
    ridgeCoef <- coef(ridgeFit, s="lambda.min")
    
    ridgeOutput[i, ] <- as.vector(ridgeCoef)
  }
  return(ridgeOutput)
}


RidgeEpan <- function(
    equation,   # The mathematical expression our data comes from
    left,       # Left endpoint of window
    right,      # Right endpoint of window
    numTerms,   # Number of Terms in the Taylor polynomial
    sigma,      # Noise for randomness in generated data
    numPoints,   # How many data points to be generated
    plot = TRUE, #If you want the data to be plotted
    seed = NULL, #Optional if you want to generate data with specific seed
    ridgeBandMulti = 4 #The factor the bandwidth returned from dpill will be multiplied by
){
  
  if(is.null(seed) == FALSE){set.seed(seed)} #Set seed if specified
  
  ##############Building Simulated Data##################
  x <- sort(runif(numPoints, min=left, max=right))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  exactY <- eval(equation)
  y <- exactY + noise
  
  #dpill returns estimated optimal bandwidth for NW and LL estimators.
  #Lasso requires a larger bandwidth and should be multiplied by a specified factor
  bandwidth <- dpill(x,y)*ridgeBandMulti 
  
  #Points at which we calculate our non-parametric estimators
  evalPoints <- seq(left, right, length.out=401)
  
  ridgeOutput <- computeRidgeOutput(x, y, bandwidth, evalPoints, numTerms)
  
  ridgeError <- findError(equation, evalPoints, ridgeOutput[,1])
  
  return(ridgeError)
}


OptimalRidge <- function(bandMultipliers = seq(1,10,1), equation, sigma, sampleSize){
  
  listOfErrors <- c()
  seed <- runif(1, min=0, max=100000)
  
  for(i in 1:length(bandMultipliers)){
    nextError <- RidgeEpan(
              equation = equation,
              left = 0,
              right = 1,
              numTerms = 10,
              sigma = sigma,
              numPoints = sampleSize,
              ridgeBandMulti = bandMultipliers[i],
              seed = seed)
    listOfErrors <- c(listOfErrors, nextError)
    print(i)
  }
  par(mfrow=c(1,1))
  plot(bandMultipliers, listOfErrors, type='l')
  return(listOfErrors)
}


ManyOptimalRidge <- function(loops, equation, sigma, sampleSize){
  
  bandMultipliers <- seq(2,10,0.5)
  
  MSEList <- vector(mode='list', length=loops)
  for(i in 1:loops){
    nextRun <- OptimalRidge(bandMultipliers, equation, sigma, sampleSize)
    MSEList[[i]] <- nextRun
  }
  par(mfrow=c(1,1))
  plot(bandMultipliers, MSEList[[1]], type="l")
  for(i in 2:loops){
    lines(bandMultipliers, MSEList[[i]])
  }
  return(MSEList)
}


#The three expressions that were tested
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))
sine <- expression(sin(5*pi*x))
bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))

Test1 <- ManyOptimalRidge(5, peak, sqrt(0.5), 250)

Test2 <- ManyOptimalRidge(5, sine, 0.5, 1000)

Test3 <- ManyOptimalRidge(5, bimodal, 0.1, 250)

Test4 <- ManyOptimalRidge(5, bimodal, 0.1, 1000)

Test1

plotBandData <- function(data, ymax){
  bandMultipliers <- seq(2,10,0.5)
  plot(bandMultipliers, data[[1]], type='l', ylim=c(0,ymax))
  for(i in 2:5){lines(bandMultipliers, data[[i]])}
}
plotBandData(Test4, 0.001)


minBandwidth <- function(data){
  listOfMins <- c()
  for(i in 1:length(data)){
    nextMinIndex <- which.min(data[[i]])
    nextMinBandwidth <- 2 + nextMinIndex * 0.5
    listOfMins <- c(listOfMins, nextMinBandwidth)
  }
  return(mean(listOfMins))
}
minBandwidth(Test1)
minBandwidth(Test2)
minBandwidth(Test3)
minBandwidth(Test4)
