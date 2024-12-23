#Function that builds simulated data
buildData <- function(n, equation, sigma){
  x <- sort(runif(n, min=0, max=1))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  trueY <- eval(equation)
  y <- trueY + noise
  return(list(x=x,y=y,trueY=trueY))
}
library(locpol)
library(KernSmooth)
#True Functions----
derivCalc <- function(func, numDeriv){
  drvList <- c(func)
  nextDrv <- func
  for(i in 1:numDeriv){
    nextDrv <- D(nextDrv, 'x')
    drvList <- c(drvList, nextDrv)
  }
  return(drvList)
}
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
derivList <- derivCalc(func = peak, numDeriv = 5)
x <- seq(0,1,length.out=401)
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),eval(derivList[6]))


runExperiment <- function(epochs, n, sigma, p, bandDerv){
  
  peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
  simulatedData <- buildData(n, peak, sigma)
  
  simData <- data.frame(matrix(numeric(401*epochs), nrow=epochs))
  
  
  for(i in 1:epochs){
    bandwidth <- thumbBw(simulatedData$x, simulatedData$y, bandDerv, EpaK)
    lassoOutput <- lassoLPR(simulatedData$x, simulatedData$y, bandwidth, p)
    
    simData[i,] <-lassoOutput$lasso[,2]
  }
  
  return(simData)
}


runLocpoly <- function(epochs, n, sigma, p, bandDerv){
  
  peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
  simulatedData <- buildData(n, peak, sigma)
  simData <- data.frame(matrix(numeric(401*epochs), nrow=epochs))
  
  for(i in 1:epochs){
    bandwidth <- thumbBw(simulatedData$x, simulatedData$y, bandDerv, gaussK)
    locOut <- locpoly(simulatedData$x, simulatedData$y, drv=p, bandwidth=bandwidth)
    
    simData[i,] <- locOut$y
  }
  return(simData)
}

#Parameters not changing between Sims
epochs <- 10
n <- 500
sigma <- sqrt(0.5)
p <- 9

#Parameters that do change
bandDerv <- 3


testLocExper <- runLocpoly(epochs, n, sigma, p=1, 1)
testLocExper2 <- runLocpoly(epochs, n, sigma, p=1, 3)
testExper <- runExperiment(epochs, n, sigma, p, bandDerv)

gridPoints <- seq(0,1, length.out=401)
plot(gridPoints, trueFunctions[[2]], type='l')
lines(gridPoints, testExper[2,], col='red')
lines(gridPoints, testLocExper[2,], col='blue')



