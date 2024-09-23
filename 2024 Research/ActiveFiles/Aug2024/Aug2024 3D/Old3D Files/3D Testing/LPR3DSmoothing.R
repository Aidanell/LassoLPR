library(KernSmooth)
library(graphics)
library(grDevices)
library(mvtnorm)
library(glmnet)
library(scatterplot3d)
library(interp)

#This file produces an iteration of the LPR Smoothing method in 3 dimensions
#It outputs results that can plot up to the 3rd derivatives of the predictons
#and compare with desired result and the locpoly method using the PlotData
#function.

#This only works for a single equation and cannot be swapped without manually
#changing equations to compare with

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


PlotData <- function(num, limits = NULL){
  #Plotting
  par(mfrow=c(2,2))
  if(is.null(limits)){
    persp(trueX, trueY, trueResult[[num]], col='lightblue',theta=30, phi=20,
          ticktype='detailed', shade=0.3, main="Exact Graph")
    persp(evalPoints, evalPoints, lassoResult[[num]], col='lightblue',theta=30, phi=20,
          ticktype='detailed', shade=0.3, main="Lasso Result")
    persp(evalPoints, evalPoints, locpolyResult[[num]], col='lightblue',theta=30, phi=20,
          ticktype='detailed', shade=0.3, main="LocPoly Result")
  }else{
    persp(trueX, trueY, trueResult[[num]], col='lightblue',theta=30, phi=20,
          ticktype='detailed', shade=0.3, main="Exact Graph", zlim=limits)
    persp(evalPoints, evalPoints, lassoResult[[num]], col='lightblue',theta=30, phi=20,
          ticktype='detailed', shade=0.3, main="Lasso Result", zlim=limits)
    persp(evalPoints, evalPoints, locpolyResult[[num]], col='lightblue',theta=30, phi=20,
          ticktype='detailed', shade=0.3, main="LocPoly Result", zlim=limits)
  }
  
}


CalculateErrors <- function(trueResult, lassoResult, locpolyResult){
  
  LassoError <- vector(length=10)
  LocpolyError <- vector(length=10)
  
  for(i in 1:10){
    LassoError[i] <- sum((lassoResult[[1]] - trueResult[[i]])**2) / length(lassoResult[[1]])
    LocpolyError[i] <- sum((locpolyResult[[1]] - trueResult[[i]])**2) / length(lassoResult[[1]])
  }
  LassoDoesBetter <- LassoError < LocpolyError
  AllErrorData <- data.frame("LassoMSE" = LassoError, "LocpolyMSE" = LocpolyError,
                             "LassoBetter" = LassoDoesBetter)
  return(AllErrorData)
}


#Functions for comparing derivatives----
f <- function(x,y){return(2-5*x-5*y +5*exp(-20*(x-0.5)**2 -20*(y-0.5)**2))}

fx <- function(x,y){return(-200*(x-0.5)*exp(-20*(x-0.5)**2 -20*(y-0.5)**2) - 5)}
fxx <- function(x,y){return((8000*(x**2) - 8000*x + 1800)*exp(-20*(x-0.5)**2 - 20*(y-0.5)**2))}
fxxx <- function(x,y){return(-1*(320000*(x**3) - 480000*(x**2) + 216000*x - 28000)*exp(-20*(x-0.5)**2 -20*(y-0.5)**2))}

fy <- function(x,y){return(-200*(y-0.5)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2) - 5)}
fyy <- function(x,y){return((8000*(y**2) - 8000*y + 1800)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2))}
fyyy <- function(x,y){return(-1*(320000*(y**3) - 480000*(y**2) + 216000*y - 28000)*exp(-20*(y-0.5)**2 -20*(x-0.5)**2))}

fxy <- function(x,y){return(2000*(2*x - 1)*(2*y - 1)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2))}
fxxy <- function(x,y){return(-4000*(40*(x**2) - 40*x +9)*(2*y - 1)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2))}
fxyy <- function(x,y){return(-4000*(40*(y**2) - 40*y +9)*(2*x - 1)*exp(-20*(y-0.5)**2 - 20*(x-0.5)**2))}

#Script Parameters----
length <- 40
xylim <- c(0,1)
numPoints <- 500
sigma <- sqrt(0.5)

#Defining matricies of true data for the function and all its derivatives----
trueX <- trueY <- seq(xylim[1], xylim[2], length.out = length)
trueResult <- vector(mode='list', length=10)
trueResult[[1]] <- outer(trueX, trueY, f)
trueResult[[2]] <- outer(trueX, trueY, fx)
trueResult[[3]] <- outer(trueX, trueY, fy)
trueResult[[4]] <- outer(trueX, trueY, fxx)
trueResult[[5]] <- outer(trueX, trueY, fxy)
trueResult[[6]] <- outer(trueX, trueY, fyy)
trueResult[[7]] <- outer(trueX, trueY, fxxx)
trueResult[[8]] <- outer(trueX, trueY, fxxy)
trueResult[[9]] <- outer(trueX, trueY, fxyy)
trueResult[[10]] <- outer(trueX, trueY, fyyy)

#Script----
#x/y values we evaluate the function at
evalPoints <- trueX

#Building Data
x <- sort(runif(numPoints, min=xylim[1], max=xylim[2]))
y <- runif(numPoints, min=xylim[1], max=xylim[2])
y <- y[order(match(y, x))]
z <- f(x, y) + rnorm(numPoints, sd=sigma)

generatedData <- data.frame("x"=x, "y"=y, "z"=z)


#Bandwidths
Xbandwidth <- dpill(x,z) * 6
Ybandwidth <- dpill(y,z) * 6

#Building Nonparametric models

lassoResult <- buildLassoMatricies(x, y, z, evalPoints, numPoints, length, Xbandwidth, Ybandwidth)
domainSize <- xylim[2] - xylim[1]
locpolyResult <- buildLocPolyMatricies(x,y,z, degree=3, Xbandwidth/domainSize, Ybandwidth/domainSize)
ErrorDf <- CalculateErrors(trueResult, lassoResult, locpolyResult)

#Creating object that holds all data from simulation
LPRResults <- vector(mode="list", 0)
LPRResults$generatedData <- generatedData
LPRResults$ErrorDf <- ErrorDf
LPRResults$lassoResult <- lassoResult
LPRResults$locpolyResult <- locpolyResult

#Defining xylim's to better see plots
flim <- c(-8,2.1)
fxlim <- c(-25, 15)
fylim <- c(-25,15)
fxxxlim <- c(-1500,1500)
PlotData(2, flim)
LPRResults


persp(evalPoints, evalPoints, lassoResult[[2]], col='lightblue',theta=30, phi=20,
      ticktype='detailed', shade=0.3, main="Lasso Result")
