

#Create neccesary data----
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2)) # sigma = sqrt(0.5)
p <- 10
sigma <- sqrt(0.5)
deriv <- 1 #Used for bandwidth calculations
sampleData <- buildData(n, line, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,EpaK)

#Function that runs lasso in a single function given data----
LassoLPR <- function(x, y, p, bandwidth){
  
  gridPoints <- seq(0, 1, length.out=401)
  lassoOutput <- matrix(nrow=401, ncol=p)
  ListOfLambdas <- c()
  
  #Estimate lasso polynomial at every gridpoint
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    X <- buildFeature(currentPoint, p, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    lassoFit <- cv.glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
    ListOfLambdas <- c(ListOfLambdas, lassoFit$lambda.min)
    
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  
  #Smooth Lambdas output with loess, and refit with new lambda graph. This reduces variance
  smoothedLambas <- lowess(gridPoints, ListOfLambdas, f=1/10)$y
  smoothLassoOutput <- matrix(nrow=401, ncol=p)
  
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    
    X <- buildFeature(currentPoint, p, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    lassoFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
    
    # Selects the desired coefficients with penalty lambda
    smoothLassoCoef <- coef(lassoFit, s=smoothedLambas[i])
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef) 
  }
  
  #Return list of smooth and unsmooth lasso, alongside the lambdas
  return(list(lasso=lassoOutput, smoothLasso=smoothLassoOutput,
              unsmoothedLambdas=ListOfLambdas, smoothedLambdas=smoothedLambas))
}


#Function that uses a specified number of lambdas----
#Should reduct time, may worsen results
NumLambdaLassoLPR <- function(x, y, p, bandwidth, numlambdas){
  
  gridPoints <- seq(0, 1, length.out=401)
  lassoOutput <- matrix(nrow=401, ncol=p)
  ListOfLambdas <- c()
  
  #Estimate lasso polynomial at every gridpoint
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    X <- buildFeature(currentPoint, p, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    lassoFit <- cv.glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7, nlambda=numlambdas)
    ListOfLambdas <- c(ListOfLambdas, lassoFit$lambda.min)
    
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  
  #Smooth Lambdas output with loess, and refit with new lambda graph. This reduces variance
  smoothedLambas <- lowess(gridPoints, ListOfLambdas, f=1/10)$y
  smoothLassoOutput <- matrix(nrow=401, ncol=p)
  
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    
    X <- buildFeature(currentPoint, p, x)
    lassoWeights <- computeWeights(x, currentPoint, bandwidth)
    lassoFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
    
    # Selects the desired coefficients with penalty lambda
    smoothLassoCoef <- coef(lassoFit, s=smoothedLambas[i])
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef) 
  }
  
  #Return list of smooth and unsmooth lasso, alongside the lambdas
  return(list(lasso=lassoOutput, smoothLasso=smoothLassoOutput,
              unsmoothedLambdas=ListOfLambdas, smoothedLambdas=smoothedLambas))
}



set.seed(202)
start1 <- Sys.time()
control <- LassoLPR(sampleData$x, sampleData$y, p, bandwidth)
end1 <- Sys.time()

set.seed(202)
start2 <- Sys.time()
lessLambdas <- NumLambdaLassoLPR(sampleData$x, sampleData$y, p, bandwidth, 10)
end2 <- Sys.time()


time1 <- end1 - start1
time2 <- end2 - start2
print(time1)
print(time2)

gridPoints <- seq(0,1, length.out=401)

plot(seq(0,1, length.out=401), control$smoothLasso[,1], type='l', col='blue')
lines(seq(0,1, length.out=401), lessLambdas$smoothLasso[,1], type='l', col='orange')


plot(seq(0,1, length.out=401), control$lasso[,1], type='l', col='blue')
lines(seq(0,1, length.out=401), lessLambdas$lasso[,1], type='l', col='orange')

plot(gridPoints, control$unsmoothedLambdas, type='l', col='blue')
lines(gridPoints, lessLambdas$unsmoothedLambdas, col='orange')

plot(gridPoints, control$smoothedLambdas, type='l', col='blue')
lines(gridPoints, lessLambdas$smoothedLambdas, col='orange')
