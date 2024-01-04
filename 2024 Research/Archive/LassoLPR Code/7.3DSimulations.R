library(KernSmooth)
library(graphics)
library(grDevices)
library(mvtnorm)
library(glmnet)
library(scatterplot3d)
library(interp)

#Helper Functions----
build2DFeature <- function(midX, midY, xValues, yValues){
  
  X <- matrix(nrow=length(xValues), ncol=9)
  
  X[,1] <- (xValues - midX)
  X[,2] <- (yValues - midY)
  X[,3] <- ((xValues - midX)**2)/2
  X[,4] <- (xValues - midX)*(yValues - midY)
  X[,5] <- ((yValues - midY)**2)/2
  X[,6] <- ((xValues - midX)**3)/6
  X[,7] <- ((xValues - midX)**2 * (yValues - midY)) /2
  X[,8] <- ((xValues - midX) * (yValues - midY)**2) /2
  X[,9] <- ((yValues - midY)**3)/6
  
  return(X)
}

epan <- function(x, bandwidth){
  return(3/4 * (1-x**2))
}

computeWeights <- function(x, midpoint, bandwidth){
  output <- vector(length=length(x))
  for(i in 1:length(x)){
    
    diff <- abs(x[i]-midpoint) / bandwidth
    if(diff > 1){output[i] <- 0}
    else{output[i] <- epan(diff)/ (bandwidth)}
    
  }
  return(output)
}

buildLasso <- function(x,y,z, evalPoints, numPoints){
  
  #Bandwidths
  Xbandwidth <- dpill(x,z) * 8
  Ybandwidth <- dpill(y,z) * 8
  
  LassoMatrix <- matrix(nrow=length, ncol=length)
  
  for(i in 1:length){
    currentX <- evalPoints[i]
    for(j in 1:length){
      currentY <- evalPoints[j]
      
      X <- build2DFeature(currentX, currentY, x, y)
      
      xWeights <- computeWeights(x, currentX, Xbandwidth)
      yWeights <- computeWeights(y, currentY, Ybandwidth)
      weights <- xWeights * yWeights
      
      lassoFit <- cv.glmnet(X, z, weights = weights, standardize=TRUE)
      lassoCoef <- as.vector(coef(lassoFit, s="lambda.min"))
      
      LassoMatrix[i,j] <- lassoCoef[1]
    }
    print(i)
  }
  return(LassoMatrix)
}
#Script----
f <- function(x,y){exp(-0.5*(x**2 + y**2))}
length <- 40
xylim <- c(-1,1)
numPoints <- 100
sigma <- 0.1

exactX <- exactY <- seq(xylim[1], xylim[2], length.out = length)
exactZ <- outer(exactX,exactY,f)
evalPoints <- exactX


#Building Data
x <- sort(runif(numPoints, min=xylim[1], max=xylim[2]))
y <- runif(numPoints, min=xylim[1], max=xylim[2])
y <- y[order(match(y, x))]
z <- f(x, y) + rnorm(numPoints, sd=sigma)


#Building Nonparametric models
lassoModel <- buildLasso(x, y, z, evalPoints, numPoints)
domainSize <- xylim[2] - xylim[1]
LLResult <- interp::locpoly(x,y,z, degree=1, h=c(Xbandwidth/domainSize, Ybandwidth/domainSize))
NWResult <- interp::locpoly(x,y,z, degree=0, h=c(Xbandwidth/domainSize, Ybandwidth/domainSize))


#Computing Errors
lassoErr <- mean((lassoModel - exactZ)**2)
NWErr <- mean((NWResult$z - exactZ)**2)
LLErr <- mean((LLResult$z - exactZ)**2)
print(lassoErr)
print(NWErr)
print(LLErr)

#Plotting
par(mfrow=c(2,2))
persp(exactX, exactY, exactZ, col='lightblue',theta=30, phi=20,
      ticktype='detailed', shade=0.3, main="Exact Graph", zlim=c(0.3,1.1))
persp(evalPoints, evalPoints, lassoModel, col='lightblue',theta=30, phi=20,
      ticktype='detailed', shade=0.3, main="Lasso Result", zlim=c(0.3,1.1))
persp(NWResult$x, NWResult$y, NWResult$z, col='lightblue',theta=30, phi=20,
      ticktype='detailed', shade=0.3, main="Nadaraya-Watson Result", zlim=c(0.3,1.1))
persp(LLResult$x, LLResult$y, LLResult$z, col='lightblue',theta=30, phi=20,
      ticktype='detailed', shade=0.3, main="Local Linear Result", zlim=c(0.3,1.1))

