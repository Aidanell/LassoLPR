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

buildLassoMatricies <- function(x,y,z, evalPoints, numPoints, length, Xbandwidth, Ybandwidth){
  
  
  ListOfLassoMatricies <- vector(mode='list', length=10)
  for(i in 1:10){ListOfLassoMatricies[[i]] <- matrix(nrow=length, ncol=length)}
  
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
      
      for(k in 1:10){
        ListOfLassoMatricies[[k]][i,j] <- lassoCoef[k]
      }
    }
    print(i)
  }
  return(ListOfLassoMatricies)
}


buildLocPolyMatricies <- function(x,y,z, degree, Xbandwidth, Ybandwidth){
  LLResult <- interp::locpoly(x,y,z, degree=degree, pd='all',
                  h=c(Xbandwidth, Ybandwidth))
  
  ListOfLocPolyMatricies <- vector(mode='list', length=10)
  ListOfLocPolyMatricies[[1]] <- LLResult$z
  ListOfLocPolyMatricies[[2]] <- LLResult$zx
  ListOfLocPolyMatricies[[3]] <- LLResult$zy
  ListOfLocPolyMatricies[[4]] <- LLResult$zxx
  ListOfLocPolyMatricies[[5]] <- LLResult$zxy
  ListOfLocPolyMatricies[[6]] <- LLResult$zyy
  ListOfLocPolyMatricies[[7]] <- LLResult$zxxx
  ListOfLocPolyMatricies[[8]] <- LLResult$zxxy
  ListOfLocPolyMatricies[[9]] <- LLResult$zxyy
  ListOfLocPolyMatricies[[10]] <- LLResult$zyyy
  return(ListOfLocPolyMatricies)
}

#Script----
f <- function(x,y){exp(-0.5*(x**2 + y**2))}
length <- 40
xylim <- c(-1,1)
numPoints <- 250
sigma <- 0.1

exactX <- exactY <- seq(xylim[1], xylim[2], length.out = length)
exactZ <- outer(exactX,exactY,f)
evalPoints <- exactX


#Building Data
x <- sort(runif(numPoints, min=xylim[1], max=xylim[2]))
y <- runif(numPoints, min=xylim[1], max=xylim[2])
y <- y[order(match(y, x))]
z <- f(x, y) + rnorm(numPoints, sd=sigma)
DataList <- vector(mode='list', length=9)

#Bandwidths
Xbandwidth <- dpill(x,z) * 8
Ybandwidth <- dpill(y,z) * 8

#Building Nonparametric models
lassoModel <- buildLassoMatricies(x, y, z, evalPoints, numPoints, length, Xbandwidth, Ybandwidth)
domainSize <- xylim[2] - xylim[1]
NWResult <- buildLocPolyMatricies(x,y,z, degree=3, Xbandwidth/domainSize, Ybandwidth/domainSize)

thirdYDeriv <- function(x,y){
  return(log(2)**2 * y * (2*log(2)*(y**2 - 3))* 2**(-y**2 -x**2 +2))
}
thirdDerivZ <- outer(exactX, exactY, thirdYDeriv)
PlotData <- function(num){
  #Plotting
  par(mfrow=c(2,2))
  persp(exactX, exactY, exactZ, col='lightblue',theta=30, phi=20,
        ticktype='detailed', shade=0.3, main="Exact Graph", zlim=c(0.3,1.1))
  persp(evalPoints, evalPoints, lassoModel[[num]], col='lightblue',theta=30, phi=20,
        ticktype='detailed', shade=0.3, main="Lasso Result", zlim=c(0.3,1.1))
  persp(evalPoints, evalPoints, NWResult[[num]], col='lightblue',theta=30, phi=20,
        ticktype='detailed', shade=0.3, main="Nadaraya-Watson Result", zlim=c(0.3,1.1))
  
}

PlotData(2)

differenceMatrixLasso <- (exactZ - lassoModel[[1]])**2
differenceMatrixLocPoly <- (exactZ - NWResult[[1]])**2
print(sum(differenceMatrixLasso))
print(sum(differenceMatrixLocPoly))






