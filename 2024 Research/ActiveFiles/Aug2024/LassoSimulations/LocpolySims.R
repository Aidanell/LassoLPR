
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2)) # sigma = sqrt(0.5)
equation <- peak
n <- 500
p <- 10
sigma <- sqrt(0.5)
gridPoints <- seq(0,1,length.out=401)

locOutput <- locpolyEpochs(500, n, equation, sigma, deriv=0)
aveLocOut <- colMeans(locOutput)

plot(gridPoints, trueFunctions[[1]])
lines(gridPoints, aveLocOut)


locpolyEpochs <- function(epochs, n, equation, sigma, deriv, degree = deriv+1){
  locOutput <- rep(0,401)
  for(i in 1:epochs){
    data <- buildData(n, equation, sigma)
    bandwidth <- thumbBw(data$x, data$y, deriv, gaussK)
    locEpoch <- locpoly(data$x, data$y, drv=deriv, degree=degree, bandwidth=bandwidth)
    locOutput <- rbind(locOutput, locEpoch$y)
  }
  locOutput <- locOutput[-1,] #Drop first row of 0's
  return(locOutput)
}

calcLocpolyVar <- function(locpolyData, trueOutputs, epochs){
  varData <- matrix(0,nrow(locpolyData), length(trueOutputs))
  
  for(i in 1:nrow(locpolyData)){
    varData[i,] <- (locpolyData[i,] - trueOutputs)**2
  }
  return(colMeans(varData))
}



testLocData <- buildData(n, equation, sigma)
testLocBand <- thumbBw(testLocData$x, testLocData$y, 0, gaussK)
testLocOut <- locpoly(testLocData$x, testLocData$y, degree=1, bandwidth=testLocBand)

plot(gridPoints, trueFunctions[[1]], type='l')
lines(testLocOut)



