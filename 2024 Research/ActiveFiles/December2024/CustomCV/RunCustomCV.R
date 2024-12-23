
library(locpol)
library(KernSmooth)


#----Pre-Setup----
#Function that builds simulated data
buildData <- function(n, equation, sigma){
  x <- sort(runif(n, min=0, max=1))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  trueY <- eval(equation)
  y <- trueY + noise
  return(list(x=x,y=y,trueY=trueY))
}

#Create neccesary data----
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5*x)
p <- 9
sigma <- sqrt(0.5)
deriv <- 1 #Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,EpaK)


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

#Builds the true datapoints for the function. Helps with MSE
derivList <- derivCalc(func = peak, numDeriv = 5)
x <- seq(0,1,length.out=401)
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),eval(derivList[6]))






testingCustomCV <- function(x, y, h, p=10, numGridPoints=401){
  
}