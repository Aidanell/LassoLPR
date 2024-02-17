library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)

plotLassoObject <- function(lassoObject, deriv = 0, lassoOnly = TRUE, points=FALSE){
  
  D0 <- lassoObject[[deriv + 1]]
  x <- D0$EvalPoints

  if(points == FALSE | deriv > 0){
    plot(D0$EvalPoints, eval(D0$DerivFunc), main=paste("Derivative",as.character(deriv)), ylab='y', type='l', lwd=3)
  }else{
    plot(D0$x, D0$y, main=paste("Derivative",as.character(deriv)), ylab='y', lwd=3, col=rgb(0,0,1,0.2), pch=16)
    lines(D0$EvalPoints, eval(D0$DerivFunc))
  }
  lines(D0$LassoFit, col='red')
  if(lassoOnly == FALSE){
    lines(D0$RidgeFit, col='green')
    lines(D0$LocpolyFit, col='blue')
  }
}


#Parameters 
poly <- expression(x**3 - 1.5*x**2 + 0.7)
left <- -0.6
right = 1.4
numTerms = 10
sigma = 0.3
numPoints = 100
seed = 2000

#Run Lasso using cross validated lambdas
UnsmoothedData <- UnsmoothedLambdasLPR(equation = poly, left=left, right=right,
                                       numTerms=numTerms, sigma=sigma,
                                       numPoints=numPoints, seed=seed)
#Performing loess to smoothen lambdas
smoothedLambas <- lowess(UnsmoothedData[[1]]$EvalPoints, UnsmoothedData[[5]], f=1/10)$y
#Applying smoothed lambdas onto method
SmoothedData <- SmoothedLambdasLPR(equation = poly, left=left, right=right,
                                   numTerms=numTerms, sigma=sigma,
                                   numPoints=numPoints, seed=seed,
                                   lambdas = smoothedLambas)



#Plotting Comparisons of Un-Smoothed vs Smoothed Lambdas
#READ FOR PLOTTING OPTIONS
#-deriv is what derivative will display (0,1,2)
#-Changing lassoOnly to False will show ridge and locpoly methods as well
#-Changing points = TRUE will display data points (must be deriv=0)
par(mfrow=c(4,2))
for(i in 0:2){
  plotLassoObject(UnsmoothedData, deriv=i, points=TRUE, lassoOnly = FALSE)
  plotLassoObject(SmoothedData, deriv=i, points=TRUE, lassoOnly = FALSE)
}
#Plotting lambda magnitudes
plot(UnsmoothedData[[1]]$EvalPoints, UnsmoothedData[[5]], main="Lambda Magnitude")
lambdaLimits <- c(min(UnsmoothedData[[5]]), max(UnsmoothedData[[5]]))
plot(SmoothedData[[1]]$EvalPoints, SmoothedData[[5]], main="Lambda Magnitude", ylim=lambdaLimits)

