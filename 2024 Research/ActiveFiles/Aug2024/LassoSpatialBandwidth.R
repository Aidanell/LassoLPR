
#This file is intended to calculate how lasso and ridge when using a similar
#rule of thumb as used in Siefert and Gasser 1998.
#Both smoothed and unsmoothed lambdas will be compared

#If lasso does better with the same lambda used in ridge rule of thumb
#then lasso should be better

numTerms <- 10
numPoints = 100
set.seed(4)
bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2)) #sigma = 0.1
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2)) # sigma = sqrt(0.5)
sine <- expression(sin(5*pi*x)) #sigma = 0.5

equation <- peak
sigma <- sqrt(0.5)

#Helper Methods----

# This function builds our feature matrix.
buildFeature <- function(midpoint, numTerms, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = numTerms)
  
  for(i in 1:numTerms){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}


#Computes MSE
findError <- function(equation, estimateX, estimateY){
  x <- estimateX
  trueY <- eval(equation)
  MSE <- mean((estimateY - trueY)^2)
  return(MSE)
}



#Excludes the edges of the function in the MSE calculation
findInteriorError <- function(equation, estimateX, estimateY, percentageCut=0.05){
  leftCutoff = as.integer(length(estimateX) * percentageCut)
  rightCutoff = as.integer(length(estimateY) - leftCutoff)
  
  x <- estimateX[leftCutoff : rightCutoff]
  y <- estimateY[leftCutoff : rightCutoff]
  trueY <- eval(equation)
  MSE <- mean((y - trueY)^2)
  return(MSE)
}



#Calculates true derivatives needed to compare estimates of derivatives with real
derivCalc <- function(func, numDeriv){
  drvList <- c(func)
  nextDrv <- func
  for(i in 1:numDeriv){
    nextDrv <- D(nextDrv, 'x')
    drvList <- c(drvList, nextDrv)
  }
  return(drvList)
}


#Function of the Epanechnikov kernel
epan <- function(x, bandwidth){
  return(3/4 * (1-x**2))
}


#Computes the scaled kernel weights for every input x
computeWeights <- function(x, midpoint, bandwidth){
  output <- vector(length=length(x))
  for(i in 1:length(x)){
    
    diff <- abs(x[i]-midpoint) / bandwidth
    if(diff > 1){output[i] <- 0}
    else{output[i] <- epan(diff)/ (bandwidth)}
    
  }
  return(output)
}


#Calculates all the derivatives needed for plotting exact functions and calculating errors
derivList <- derivCalc(func = equation, numDeriv = numTerms)

##############Building Simulated Data##################


x <- sort(runif(numPoints, min=left, max=right))
noise <- rnorm(n = length(x), mean = 0, sd = sigma)
exactY <- eval(equation)
y <- exactY + noise


####Bandwidth Calculations----
#For Epanechnikov kernel, max(K(u)) = 3/4 and R(K) = 1/2
#Therefore bandwidth = 5*delta)/16

locpolyBandwidth <- dpill(x,y)
lassoBandwidth <- locpolyBandwidth*8
evalPoints <- seq(left, right, length.out=401)


##Computing Lasso----

lassoOutput <- matrix(nrow=401, ncol=numTerms+1)
ridgeOutput <- matrix(nrow=401, ncol=numTerms+1)

lambdaTracking <- c()
XTildeTracking <- c()

for(i in 1:length(evalPoints)){
  currentPoint <- evalPoints[i]
  
  X <- buildFeature(currentPoint, numTerms, x)
  
  
  lassoWeights <- computeWeights(x, currentPoint, lassoBandwidth)
  XTilde <- sum(lassoWeights * x) / sum(lassoWeights)
  
  lassoFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
  ridgeFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=0, maxit=10**7)
  
  RidgeThumb = 5*abs(currentPoint - XTilde)/(16*lassoBandwidth)
  
  # Selects the desired coefficients with penalty lambda
  lassoCoef <- coef(lassoFit, s=RidgeThumb)
  ridgeCoef <- coef(ridgeFit, s=RidgeThumb)
  
  
  #The i'th row contains all derivative estimation at that specificed evalpoint
  lassoOutput[i, ] <- as.vector(lassoCoef)
  ridgeOutput[i, ] <- as.vector(ridgeCoef) 
  
  lambdaTracking <- c(lambdaTracking, RidgeThumb)
  XTildeTracking <- c(XTildeTracking, XTilde)
}



#####Spatially Adaptive Ridge
#Calculate optimal ridge parameter at each point then refit with lasso

SpatialRidge <- GetSpatialRidgeList(x, evalPoints, bandwidth, XTildeTracking, ridgeOutput[,2])

SpatialLassoOutput <- matrix(nrow=401, ncol=numTerms+1)
SpatialRidgeOutput <- matrix(nrow=401, ncol=numTerms+1)

#Now reevaluate lasso and ridge with new lambdas.
for(i in 1:length(evalPoints)){
  
  currentPoint <- evalPoints[i]
  
  X <- buildFeature(currentPoint, numTerms, x)
  
  
  lassoWeights <- computeWeights(x, currentPoint, lassoBandwidth)
  XTilde <- sum(lassoWeights * x) / sum(lassoWeights)
  
  SpatialLassoFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
  SpatialRidgeFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=0, maxit=10**7)
  
  SpatialLambda = SpatialRidge[i]
  
  # Selects the desired coefficients with penalty lambda
  SpatialLassoCoef <- coef(SpatialLassoFit, s= SpatialLambda)
  SpatialRidgeCoef <- coef(SpatialRidgeFit, s= SpatialLambda)
  
  
  #The i'th row contains all derivative estimation at that specificed evalpoint
  SpatialLassoOutput[i, ] <- as.vector(SpatialLassoCoef)
  SpatialRidgeOutput[i, ] <- as.vector(SpatialRidgeCoef) 
}



plot(x,y, main="Lasso (Red) vs Ridge (Blue) Thumb")
lines(evalPoints, lassoOutput[,1], col='red', lwd=2)
lines(evalPoints, ridgeOutput[,1], col='blue', lwd=2)
lines(evalPoints, testLassoEpan$smoothLassoOutput[,1], col="orange", lwd=2)
lines(evalPoints, testLassoEpan$lassoOutput[,1], col="purple", lwd=2)

plot(x,y, main='Lasso Spatial (Red) vs Lasso Thumb (Blue)')
lines(evalPoints, SpatialLassoOutput[,1], col='red', lwd=3)
lines(evalPoints, lassoOutput[,1], col='blue', lwd=3)

plot(x,y, main='Ridge Spatial (Red) vs Ridge Thumb (Blue)')
lines(evalPoints, SpatialRidgeOutput[,1], col='red', lwd=3)
lines(evalPoints, ridgeOutput[,1], col='blue', lwd=3)

plot(evalPoints, SpatialRidge, ylim=c(0,3), main="Spatial Lambdas (Red) vs Thumb Lambdas (Blue)", col='red')
plot(evalPoints, lambdaTracking, main="Thumb Lambdas", col="blue")
points(evalPoints, testLassoEpan$lambdas$smoothed, col='orange')
points(evalPoints, testLassoEpan$lambdas$unsmoothed, col='purple')



#temp <- x
#x <- evalPoints
#deriv <- 2
#plot(evalPoints, eval(derivList[deriv+1]), main="Lasso Spatial (Red) vs Ridge Spatial (Blue), 1st Derivative")
#lines(evalPoints, SpatialLassoOutput[,deriv+1], col='red')
#lines(evalPoints, SpatialRidgeOutput[,deriv+1], col='blue')

#plot(evalPoints, eval(derivList[deriv+1]), main="Lasso Thumb (Red) vs Ridge Thumb (Blue), 1st Derivative")
#lines(evalPoints, lassoOutput[,deriv+1], col='red')
#lines(evalPoints, ridgeOutput[,deriv+1], col='blue')
#x <- temp







