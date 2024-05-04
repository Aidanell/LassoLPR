library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)



#Computes MSE
findError <- function(equation, estimateX, estimateY){
  x <- estimateX
  trueY <- eval(equation)
  MSE <- mean((estimateY - trueY)^2)
  return(MSE)
}


#The three expressions that were tested
sine <- expression(sin(5*pi*x))
bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))

#Parameters
numPoints <- 250
left <- 0
right <- 1
sigma <- sqrt(0.5)
evalPoints <- seq(left, right, length.out=401)

x <- sort(runif(numPoints, min=left, max=right))
noise <- rnorm(n = length(x), mean = 0, sd = sigma)
exactY <- eval(peak)
y <- exactY + noise

bandwidth <- dpill(x,y)
bandwidthMultipliers <- seq(0.25,15,0.25) #Bandwidths to evaluate LPR
MSEOutput <- c()

for(i in 1:length(bandwidthMultipliers)){
  locpolyFit <- KernSmooth::locpoly(x, y, degree = 3, drv = 0,
                                    bandwidth = bandwidth * bandwidthMultipliers[i], gridsize = 401) 
  locpolyError <- findError(peak, locpolyFit$x, locpolyFit$y)
  MSEOutput <- c(MSEOutput, locpolyError)
}


plot(bandwidthMultipliers, MSEOutput, main="LocPoly")





