
library(locpol)
library(KernSmooth)
library(zoo)

#Create neccesary data----
set.seed(100)
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5 + sin(5*x)/3)
#peak <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))
p <- 9
sigma <- sqrt(0.5)
#sigma <- 0.1
deriv <- 1#Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,gaussK)


#Builds the true datapoints for the function. Helps with MSE----
derivList <- derivCalc(func = peak, numDeriv = 10)
x <- seq(0,1,length.out=401)
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),
                      eval(derivList[6]), eval(derivList[7]),eval(derivList[8]),eval(derivList[9]),eval(derivList[10]))
#Calculate true Y for simulated x's
x <- sampleData$x
trueY <- eval(peak)


#This script focuses on minimizing the number of glmnet calls
#Evaluate lambda at using 25 different windows that do not overlap (or maybe do slightly)
#Calculate lambda using LOOCV. Then, interpolate the lambdas and calculate at each gridpoint
#with entire dataset


#When d=25, each interval is 0.04 (roughly 20 datapoints per window)
d <- 40

getPointsInWindow <- function(x, gridPoint, band){
  currentMax = gridPoint + band
  currentMin = gridPoint - band
  
  inWindow <- which(sampleData$x >= currentMin & sampleData$x <= currentMax)
  return(inWindow)
}

gridPoints <- (1:d / d) - (1/d)/2 #Where all the lambdas are being evaluated
lambdas <- numeric(length=d)
for(i in 1:d){
  
  currentGridPoint <- gridPoints[i]
  IndiciesInWindow <- getPointsInWindow(sampleData$x, currentGridPoint, 0.1)
  #print(IndiciesInWindow)
  currentY <- sampleData$y[IndiciesInWindow]
  currentX <- sampleData$x[IndiciesInWindow]
  
  
  X <- buildFeature(currentGridPoint, p, currentX)
  lassoWeights <- computeWeights(currentX, currentGridPoint, 0.05, kernel='norm')
  lassoFit <- glmnet::cv.glmnet(X, currentY, weights = lassoWeights, maxit=10**7, nfolds=length(currentY), grouped=FALSE, nlambda=20)
  
  lambdas[i] <- lassoFit$lambda.min
  
}


loessdf <- data.frame('gridPoints'=gridPoints, 'lambdas'=lambdas)
loessSmoothModel <- stats::loess(lambdas ~ gridPoints, data=loessdf, span=0.35, na.action=na.omit)
plot(lambdas)
lines(loessSmoothModel$fitted, type='l')

gridPoints <- seq(0,1,length.out=401)
smoothLambdas <- stats::predict(loessSmoothModel, newdata=gridPoints)

filledLambdas <- na.locf(smoothLambdas, fromLast = TRUE, na.rm=FALSE)
filledLambdas <- na.locf(filledLambdas)    
plot(gridPoints, smoothLambdas, type='l')

filledLambdas[filledLambdas < 0] <- 0

smoothLassoOutput <- matrix(nrow=401, ncol=p+1)
listOfZeroDerivatives <- numeric(401)

for(i in 1:length(gridPoints)){
  currentPoint <- gridPoints[i]
  
  X <- buildFeature(currentPoint, p, sampleData$x)
  lassoWeights <- computeWeights(sampleData$x, currentPoint, bandwidth, kernel = "norm")
  lassoFit <- glmnet::glmnet(X, sampleData$y, weights = lassoWeights, maxit=10**7)
  
  # Selects the desired coefficients with the best lambda
  smoothLassoCoef <- coef(lassoFit, s=filledLambdas[i], exact=TRUE, x=X, y=sampleData$y, weights=lassoWeights)
  #The i'th row contains all estimated p+1 derivatives estimation at the i'th gridpoint
  smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef)
  listOfZeroDerivatives[i] = sum(smoothLassoOutput[i,] == 0)
}


plot(gridPoints, smoothLambdas, ylim=c(0,0.2), col='red', type='l', lwd=2)
points((1:d / d) - (1/d)/2, lambdas, col='red', pch=16)
lines(gridPoints, customCVLasso[[4]], lwd=2)
points(gridPoints, customCVLasso[[2]], pch=16)

plotDeriv <- 0
plot(gridPoints, trueFunctions[[1+plotDeriv]], type='l')
lines(gridPoints,smoothLassoOutput[,1+plotDeriv], col='red')
lines(gridPoints, customCVLasso[[3]][,1+plotDeriv], col='blue')


plot(gridPoints, smoothLassoOutput[,3] != 0, type='l')
discreteVariableInclusion <- matrix(nrow=401, ncol=p+1)
for(i in 1:(p+1)){
  discreteVariableInclusion[,i] <- round(lowess(gridPoints, smoothLassoOutput[,i] != 0, f=0.15)$y)
}



####Now lets use the smoothLassoOutput to see which derivatives are zero, and leave them out
#of the locpoly estimation. Use a new build matrix to faciliate

buildFeatureExclude <- function(midpoint, p, xValues, termsToExclude){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  termsToExclude <- termsToExclude - 1
  
  #numCol <- p-length(termsToExclude)
  #currCol <- 1
  X <- matrix(nrow = length(xValues), ncol = p)
  for(i in 1:p){
    if(!(i %in% termsToExclude)){
      X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
      #currCol <- currCol + 1
    }else{
      X[,i] <- rep(0,length(xValues))
    }
  }
  return(X)
}


smoothLassoOutputExclude <- matrix(nrow=401, ncol=p+1)

for(i in 1:length(gridPoints)){
  currentPoint <- gridPoints[i]
  
  termsToExclude <- which(discreteVariableInclusion[i,] == 0)
  X <- buildFeatureExclude(currentPoint, p, sampleData$x, termsToExclude)

  lassoWeights <- computeWeights(sampleData$x, currentPoint, bandwidth, kernel = "norm")
  
  test <- lm(sampleData$y ~ X, weights=lassoWeights)
  lmCoefs <- as.vector(test$coefficients)
  lmCoefs[is.na(lmCoefs )] <- 0
  smoothLassoOutputExclude[i, ] <- lmCoefs
}


plotDeriv <- 3
plot(gridPoints, trueFunctions[[1+plotDeriv]], type='l')
lines(gridPoints,smoothLassoOutput[,1+plotDeriv], col='red')
lines(gridPoints,smoothLassoOutputExclude[,1+plotDeriv], col='purple')

#ts.plot(discreteVariableInclusion[,1])
