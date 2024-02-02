library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)

lassoLPR <- function(x, y, bandwidth, evaluatePoints, numDerivatives = 10){
  
  lassoOutput <- matrix(nrow=401, ncol=numDerivatives+1) #Plus one for constant column
  ListOfLambdas <- c()
  for(i in 1:length(evaluatePoints)){
    point <- evaluatePoints[i]
    X <- buildFeatureMatrix(point, numTerms = numDerivatives, x)
    weights <- computeWeights(x, midpoint = point, bandwidth)
    
    lassoFit <- cv.glmnet(X, y, weights = weights, standardize=TRUE, alpha=1)
    # Selects the desired coefficients with penalty lambda
    lassoCoef <- coef(lassoFit, s="lambda.min")
    ListOfLambdas <- c(ListOfLambdas, lassoFit$lambda.min)
    
    #The i'th row contains all derivative estimation at that specificed evalpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
    
  }
  
  lassoData <- vector(mode='list', 2)
  lassoData[[1]] <- lassoOutput
  lassoData[[2]] <- ListOfLambdas
  return(lassoData)
} 


#Computes the scaled kernel weights for every input x
computeWeights <- function(x, midpoint, bandwidth){
  
  #Function of the Epanechnikov kernel
  epan <- function(x, bandwidth){
    return(3/4 * (1-x**2))
  }
  
  weights <- vector(length=length(x))
  for(i in 1:length(x)){
    
    diff <- abs(x[i]-midpoint) / bandwidth
    if(diff > 1){weights[i] <- 0}
    else{weights[i] <- epan(diff)/ (bandwidth)}
  }
  
  return(weights)
}



buildFeatureMatrix <- function(midpoint, numTerms, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = numTerms)
  
  for(i in 1:numTerms){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}
#Script----
#Reading Data In
motorcycle <- data.frame(read.csv("motorcycle.txt", head = TRUE, sep=" "))

#Changing splices changes where the x-value cutoffs are for each grouping
firstSplice <- 15
secondSplice <- 28
firstSpliceIndicies <- which(motorcycle[,1] < firstSplice)
secondSpliceIndicies <- which(motorcycle[,1] < secondSplice)
firstSpliceIndex <- tail(firstSpliceIndicies, 1) #works since this column is ordered already
secondSpliceIndex <- tail(secondSpliceIndicies, 1)

#Split into groups, divide each by their standard deviation, and recombine
group1 <- motorcycle[1:firstSpliceIndex,]
group2 <- motorcycle[firstSpliceIndex:secondSpliceIndex,]
group3 <- motorcycle[secondSpliceIndex:132,]

sdg1 <- sd(group1[,2])
sdg2 <- sd(group2[,2])
sdg3 <- sd(group3[,2])

scaledData <- motorcycle
scaledData[1:firstSpliceIndex, 2] <- group1[,2] / sdg1
scaledData[firstSpliceIndex:secondSpliceIndex, 2] <- group2[,2] / sdg2
scaledData[secondSpliceIndex:132, 2] <- group3[,2] / sdg3

#Calculating bandwidth and gridpoints
bandwidth <- dpill(scaledData[,1], scaledData[,2])
lassoBandwidth <- bandwidth * 8
pointstoEvaluateAt <- seq(min(scaledData[,1]), max(scaledData[,1]), length.out=401)

#Perform Lasso
motorLasso <- lassoLPR(scaledData[,1], scaledData[,2], lassoBandwidth, pointstoEvaluateAt)

#Unstandardize lasso output by grouping
convertedData <- rep(0, 401)
for(i in 1:401){
  if(pointstoEvaluateAt[i] <= firstSplice){
    convertedData[i] <- motorLasso[[1]][i,1] * sdg1
  }else if(pointstoEvaluateAt[i] <= secondSplice){
    convertedData[i] <- motorLasso[[1]][i,1] * sdg2
  }else{
    convertedData[i] <- motorLasso[[1]][i,1] * sdg3
  }
}

#Plotting
par(mfrow=c(1,2))
plot(motorcycle, main='Lasso vs Locpoly, Scaling first')
lines(pointstoEvaluateAt, convertedData, col='red')

#Original Lasso For Comparison
bandwidth <- dpill(motorcycle[,1], motorcycle[,2])
lassoBandwidth <- bandwidth * 8
lassoBandwidth
pointstoEvaluateAt <- seq(min(motorcycle[,1]), max(motorcycle[,1]), length.out=401)
motorLassoOrig <- lassoLPR(motorcycle[,1], motorcycle[,2], lassoBandwidth, pointstoEvaluateAt)

plot(motorcycle, main='Lasso vs Locpoly')
lines(pointstoEvaluateAt, motorLassoOrig[[1]][,1], col='red')



