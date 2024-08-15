library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)

equation = expression(2-5*x +5*exp(-400*(x-0.5)**2))
left = 0
right = 1
numTerms = 10
sigma = sqrt(0.5)
numPoints = 500
seed=50

x <- sort(runif(numPoints, min=left, max=right))
noise <- rnorm(n = length(x), mean = 0, sd = sigma)
exactY <- eval(equation)
y <- exactY + noise

#dpill returns estimated optimal bandwidth for locpoly estimators.
#Lasso requires a larger bandwidth and should be multiplied by a specified factor
locpolyBandwidth <- dpill(x,y)
lassoBandMulti = 8
lassoBandwidth <- locpolyBandwidth*lassoBandMulti


#Points at which we calculate our non-parametric estimators
evalPoints <- seq(left, right, length.out=401)

##############Computing Smoothing Function#############

currentPoint <- evalPoints[48]
X <- buildFeature(currentPoint, numTerms, x)

lassoWeights <- computeWeights(x, currentPoint, lassoBandwidth)

lassoOutput <- matrix(nrow=100, ncol=numTerms+1)
for(i in 1:100){
  lassoFit <- cv.glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
  lassoCoef <- coef(lassoFit, s="lambda.min")
  lassoOutput[i, ] <- as.vector(lassoCoef) 
}



hist(lassoOutput[,5], breaks=20, freq=FALSE)
lassoOutput[,6]
