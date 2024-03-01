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
sigma = sqrt(0.5)
numPoints = 500
numExperiments <- 500

listOfExperiments <- vector(mode='list', length=numExperiments)


for(i in 1:numExperiments){
  listOfExperiments[[i]] <- RunLPRSimulation(equation = poly, left=left, right=right,
                                             numTerms=numTerms, sigma=sigma,
                                             numPoints=numPoints)
  par(mfrow=c(2,2))
  for(j in 0:2){
    plotData(listOfExperiments[[i]], deriv=j)
  }
  plot(listOfExperiments[[i]]$evalPoints, listOfExperiments[[i]]$lambdas$unsmoothed, main="Lambdas", type='l', col='red')
  lines(listOfExperiments[[i]]$evalPoints, listOfExperiments[[i]]$lambdas$smoothed, col='orange')
  
  print(i)
}


plotBiasAndVar(listOfExperiments, deriv=3)
