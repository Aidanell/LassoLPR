#This file is exploring the basics of creating new synthetic data using 
#synthetic discrete derivatives

set.seed(20)
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5 + sin(5*x)/3)
#peak <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))
p <- 9
sigma <- sqrt(0.5)
#sigma <- 0.1
deriv <- 3#Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,gaussK)

derivedData <- numeric(500)
for(i in 2:499){
  derivedData[i] <- (sampleData$y[i+1] - sampleData$y[i]) / (sampleData$x[i+1] - sampleData$x[i]) 
}

plot(sampleData$x, derivedData, ylim=c(-10**3, 10**3))
lines(gridPoints, trueFunctions[[2]])     


####Use linear combination of weights
library("doBy")


calcWeight <- function(j){
  w <- (6 * j**2) / (k*(k+1)*(2*k+1))
  return(w)
}

k <- 20
synData <- c()
for(i in (k+1):(500-k)){
  
  chosenIndicies <- which.minn(abs(sampleData$x - sampleData$x[i]),2*k)
  
  belowXIn <- chosenIndicies[1:k]
  aboveXIn <- chosenIndicies[(k+1):(2*k)]
  
  belowX <- sampleData$x[belowXIn]
  aboveX <- sampleData$x[aboveXIn]
  
  belowY <- sampleData$y[belowXIn]
  aboveY <- sampleData$y[aboveXIn]
  
  output <- 0
  for(j in 1:k){
    coef <- (aboveY[j] - belowY[j])/ (aboveX[j] - belowX[j])
    output <- output + calcWeight(j) * coef
  }
  synData <- c(synData, output)
}


plot(sampleData$x[(k+1):(500-k)], synData, col='red')
lines(gridPoints, trueFunctions[[2]])     




