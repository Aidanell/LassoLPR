#Master Lasso File
library(matlib)
library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)
#Parameters
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))
poly <- expression(80*(x-0.5)**3 - 5*(x-0.5) + 0.5)
equation = poly
left = 0
right = 1
degree = 10
sigma = 2
n = 500
numEvalPoints <- 401
seed=42
set.seed(seed)

#Build Data
observedX <- sort(runif(n, min=left, max=right))
noise <- rnorm(n = length(observedX), mean = 0, sd = sigma)

#Adding noise to true data
x <- observedX
trueY <- eval(equation)
x <- NULL
observedY <- trueY + noise

#Lasso Bandwidth
bandwidth <- dpill(observedX,observedY)*4
print("Bandwidth:")
print(bandwidth)

#Points where we estimate the function
evalPoints <- seq(left, right, length.out=numEvalPoints)


linearData <- computeLinearModel(observedX, observedY, bandwidth, evalPoints, degree)







