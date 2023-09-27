library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)


computeLassoOutput <- function(x, y, bandwidth, evalPoints, numTerms){
  
  lassoOutput <- matrix(nrow=401, ncol=numTerms+1)
  for(i in 1:length(evalPoints)){
    currentPoint <- evalPoints[i]
    
    X <- buildFeature(currentPoint, numTerms, x)
    weights <- computeWeights(x, currentPoint, bandwidth)
    
    lassoFit <- cv.glmnet(X, y, weights = weights, standardize=TRUE, alpha=1)
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  return(lassoOutput)
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


useLasso <- function(x, y){
  
  evalPoints <- seq(min(x), max(x), length.out=401)
  bandwidth <- dpill(x,y) * 8
  numTerms <- 10
  
  lassoOutput <- computeLassoOutput(x, y, bandwidth, evalPoints, numTerms)
  
  
  output <- data.frame(x=evalPoints, y=lassoOutput[,1])
  return(output)
}


df <- DAAG::bomregions2021
x <- df$Year[11:length(df$Year)]
y<- df$northAVt[11:length(df$northAVt)]
par(mfrow=c(2,2))


lassoFit <- useLasso(x,y)
plot(lassoFit$x, lassoFit$y, type='l', main="Lasso")
points(x,y)


locFit <- locpoly(x,y, bandwidth=dpill(x,y))
plot(locFit, main="LocPoly")
points(x,y)


fit <- lm(y ~ x)
plot(x,y, main="Linear Regression")
abline(fit)


loessFit <- loess(y ~ x, span=0.7)
lo.x <- seq(min(x), max(x), length.out=50)
lo.y <- predict(loessFit, lo.x)
plot(x,y, main='Loess')
lines(lo.x, lo.y)
