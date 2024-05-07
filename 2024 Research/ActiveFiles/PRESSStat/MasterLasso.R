#Master Lasso File
library(matlib)
library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)
#Parameters
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))
poly <- expression(80*(x-0.5)**3 - 5*(x-0.5) + 0.5)
equation = peak
left = 0
right = 1
degree = 10
sigma = sqrt(0.5)
n = 500
numEvalPoints <- 401
seed=65
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
bandwidth <- dpill(observedX,observedY)*16
print("Bandwidth:")
print(bandwidth)

#Points where we estimate the function
evalPoints <- seq(left, right, length.out=numEvalPoints)

#Lasso LPR
LassoPRESSData <- computePRESSLasso(observedX, observedY, bandwidth, evalPoints, degree)
lassoOutput <- LassoPRESSData[1]
ListOfLambdas <- LassoPRESSData$lambdas



#Plotting 0th Derivative
par(mfrow=c(1,1))
plot(observedX, observedY, main="Regression of 0th Derivative")
lines(observedX, trueY, col='green', lwd=2)
lines(evalPoints, LassoPRESSData$PRESSDerivatives[,1], col='blue', lwd=2)
lines(evalPoints, testLassoEpan$lassoOutput[,1], col='red', lwd=2)
lines(evalPoints, testLassoEpan$smoothLassoOutput[,1], col='orange', lwd=2)
legend(x='bottomleft', legend=c("True Data", "PRESSLASSO", "CVLasso", "SmoothedLambdaCVLasso"),
       col=c("green", "blue","red","orange"), lwd=2)

#Plotting lambdas
plot(evalPoints, testLassoEpan$lambdas$unsmoothed, col='red', main="Plot of Lasso Lambdas")
points(evalPoints, ListOfLambdas, col="blue")
points(evalPoints, testLassoEpan$lambdas$smoothed, col='orange')
legend(x='topleft', legend=c("PRESSLASSO", "CVLasso", "SmoothedLambdaCVLasso"),
       col=c("blue","red","orange"), lwd=2)

#Getting derivative of true function
Derivatives <- derivCalc(equation, degree)

#J'th derivative
j <-3
#Getting true curve
x <- evalPoints
trueFunction <- Derivatives[j+1]
y <- eval(trueFunction)
if(length(y) == 1){ #For when the true function is constant
  y <- rep(y, 401)
}

#plotting
plot(x, y, main=paste("Regression of Derivative ", j), type='l', lwd=2, col='black')
x <- NULL
lines(evalPoints, LassoPRESSData$PRESSDerivatives[,j+1], col='blue', lwd=2)
lines(evalPoints, testLassoEpan$lassoOutput[,j+1], col='red', lwd=2)
lines(evalPoints, testLassoEpan$smoothLassoOutput[,j+1], col='orange', lwd=2)
legend(x='topright', legend=c("True Data", "PRESSLASSO", "CVLasso", "SmoothedLambdaCVLasso"),
       col=c("black", "blue","red","orange"), lwd=2)

