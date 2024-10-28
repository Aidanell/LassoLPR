#Function that builds simulated data
buildData <- function(n, equation, sigma){
  x <- sort(runif(n, min=0, max=1))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  trueY <- eval(equation)
  y <- trueY + noise
  return(list(x=x,y=y,trueY=trueY))
}
library(locpol)
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
derivList <- derivCalc(func = peak, numDeriv = 5)
x <- seq(0,1,length.out=401)
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),eval(derivList[6]))


#Create neccesary data----
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5*x)
p <- 9
sigma <- sqrt(0.5)
deriv <- 1 #Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv, gaussK)
bandwidth2 <- thumbBw(sampleData$x,sampleData$y,deriv, EpaK)

set.seed(40)
testLassoLib <- lassoLPR(sampleData$x, sampleData$y, bandwidth*3, p, kernel="normal")

set.seed(40)
testLassoLib2 <- lassoLPR(sampleData$x, sampleData$y, bandwidth2, p, kernel="epak")

gridPoints <- seq(0,1,length.out=401)
#plot(gridPoints, testLassoLib$lasso[,1], type='l')
locBand <- dpill(sampleData$x, sampleData$y)
locOut <- locpoly(sampleData$x, sampleData$y, bandwidth=locBand*1.3, drv=1)
plot(gridPoints, testLassoLib2$lasso[,2], type='l', col='blue')
lines(locOut)

plot(gridPoints, testLassoLib2$lasso[,2], type='l', col='blue')
multipliers <- seq(0.5, 5, length.out=20)
locMSE <- numeric(20)
for(i in 1:length(multipliers)){
  locOut <- locpoly(sampleData$x, sampleData$y, bandwidth=locBand*multipliers[i], drv=1, degree=3)
  lines(locOut)
  locMSE[i] <- mean(abs(trueFunctions[[2]] - locOut$y)**2 )
}

plot(multipliers, locMSE)
plot(gridPoints, testLassoLib2$lasso[,2], type='l', col='blue')
locOut <- locpoly(sampleData$x, sampleData$y, bandwidth=locBand*3, drv=1, degree=3)
lines(locOut)
lines(gridPoints, trueFunctions[[2]], col='purple')

