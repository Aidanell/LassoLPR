
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2)) # sigma = sqrt(0.5)
equation <- peak
n <- 500
p <- 10
sigma <- sqrt(0.5)
gridPoints <- seq(0,1,length.out=401)

locOutput <- locpolyEpochs(500, n, equation, sigma)
aveLocOut <- colMeans(locOutput)

varLocOut <- 



locpolyEpochs <- function(epochs, n, equation, sigma){
  locOutput <- rep(0,401)
  for(i in 1:epochs){
    data <- buildData(n, equation, sigma)
    bandwidth <- thumbBw(data$x,data$y,0,gaussK)
    locEpoch <- locpoly(data$x, data$y, degree=1, bandwidth=bandwidth)
    locOutput <- rbind(locOutput, locEpoch$y)
  }
  
  locOutput <- locOutput[-1,] #Drop first row of 0's
  return(locOutput)
}

