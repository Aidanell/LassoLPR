library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)

#Parameters 
#poly <- expression(x**3 - 1.5*x**2 + 0.7)
#left <- -0.6
#right = 1.4
poly <- expression(2-5*x +5*exp(-400*(x-0.5)**2))
left <- 0
right = 1
numTerms = 10
sigma = 0.3
numPoints = 500
seed = 2000
numExperiments <- 5

listOfExperiments <- vector(mode='list', length=numExperiments)
listOfSmoothedExperiments <- vector(mode='list', length=numExperiments)

for(i in 1:numExperiments){

  
  #Run Lasso using cross validated lambdas
  listOfExperiments[[i]] <- UnsmoothedLambdasLPR(equation = poly, left=left, right=right,
                                         numTerms=numTerms, sigma=sigma,
                                         numPoints=numPoints)
  
  #Performing loess to smoothen lambdas
  smoothedLambas <- lowess(listOfExperiments[[i]][[1]]$EvalPoints, listOfExperiments[[i]][[5]], f=1/10)$y
  #Applying smoothed lambdas onto method
  listOfSmoothedExperiments[[i]] <- SmoothedLambdasLPR(equation = poly, left=left, right=right,
                                     numTerms=numTerms, sigma=sigma,
                                     numPoints=numPoints, seed=listOfExperiments[[i]][[6]],
                                     lambdas = smoothedLambas)
  
  
  par(mfrow=c(4,2))
  for(j in 0:2){
    plotLassoObject(listOfExperiments[[i]], deriv=j,  lassoOnly = FALSE)
    plotLassoObject(listOfSmoothedExperiments[[i]], deriv=j, lassoOnly = FALSE)
  }
  #Plotting lambda magnitudes
  plot(listOfExperiments[[i]][[1]]$EvalPoints, listOfExperiments[[i]][[5]], main="Lambda Magnitude")
  lambdaLimits <- c(min(listOfExperiments[[i]][[5]]), max(listOfExperiments[[i]][[5]]))
  plot(listOfSmoothedExperiments[[i]][[1]]$EvalPoints, listOfSmoothedExperiments[[i]][[5]], main="Lambda Magnitude", ylim=lambdaLimits)
  
  print(i)

}


plotBiasAndVar(listOfSmoothedExperiments, listOfExperiments, deriv=2)
