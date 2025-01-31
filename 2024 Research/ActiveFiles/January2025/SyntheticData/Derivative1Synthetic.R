#This file is experimenting with appending synthetic data and using for glmnet



#Weight suggested by paper. There's probably better solution for our use case
calcWeight <- function(j, k){
  w <- (6 * j**2) / (k*(k+1)*(2*k+1))
  return(w)
}
createSyntheticData <- function(x, y, k){
  
  
  totalData <- length(x)
  synData <- c()
  for(i in (k+1):(totalData-k)){
    
    chosenIndicies <- which.minn(abs(x - x[i]),2*k)
    
    belowXIn <- chosenIndicies[1:k]
    aboveXIn <- chosenIndicies[(k+1):(2*k)]
    
    belowX <- x[belowXIn]
    aboveX <- x[aboveXIn]
    
    belowY <- y[belowXIn]
    aboveY <- y[aboveXIn]
    
    output <- 0
    for(j in 1:k){
      coef <- (aboveY[j] - belowY[j])/ (aboveX[j] - belowX[j])
      output <- output + calcWeight(j, k) * coef
    }
    synData <- c(synData, output)
  }
  return(list(y=synData, x=x[(k+1):(totalData-k)]))
}

createSyntheticDataV2 <- function(x, y, k){
  
  calcWeightV2 <- function(belowX, aboveX){
    
    aboveX <- rev(aboveX)
    w <- (belowX - aboveX)^2
    w <- w/sum(w)
    return(w)
  }
  
  totalData <- length(x)
  synData <- c()
  for(i in (k+1):(totalData-k)){
    
    # Step 1: Get indices of values <= x
    filtered_indices <- which(x < x[i])
    # Step 2: Sort indices based on the absolute difference of values from x
    sorted_indices <- filtered_indices[order(abs(x[filtered_indices] - x[i]))]
    # Step 3: Select the first 5 indices
    belowIn <- head(sorted_indices, k)
    # Step 4: Get the corresponding values
    belowX <- x[belowIn]
    
    # Step 1: Get indices of values <= x
    filtered_indices <- which(x > x[i])
    # Step 2: Sort indices based on the absolute difference of values from x
    sorted_indices <- filtered_indices[order(abs(x[filtered_indices] - x[i]))]
    # Step 3: Select the first 5 indices
    aboveIn <- head(sorted_indices, k)
    # Step 4: Get the corresponding values
    aboveX <- x[aboveIn]
    
    aboveY <- y[aboveIn]
    belowY <- y[belowIn]
    
    coef <- (aboveY - belowY)/ (aboveX - belowX)
    weights <- calcWeightV2(belowX, aboveX)
    #print(length(coef))
    print(weights)
    output <- sum(weights * coef)
    #print(coef)
    synData <- c(synData, output)
  }
  return(list(y=synData, x=x[(k+1):(totalData-k)]))
}

set.seed(200)
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5 + sin(5*x)/3)
#peak <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))
p <- 9
sigma <- sqrt(0.5)
#sigma <- 0.1
deriv <- 3#Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,EpaK)
gridPoints <- seq(0,1, length.out=401)
numGridPoints <- length(gridPoints)

#Builds the true datapoints for the function. Helps with MSE----
derivList <- derivCalc(func = peak, numDeriv = 5)
x <- seq(0,1,length.out=401)
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),eval(derivList[6]))
#Calculate true Y for simulated x's
x <- sampleData$x
trueY <- eval(peak)


addSynthticData <- function(X, synX, j, midpoint){
  
  newX <- matrix(nrow=length(synX), ncol=p+1)
  for(i in 1:(j-1)){
    newX[,i] <- rep(0, length(synX))
  }
  newX[,j+1] <- rep(1, length(synX))
  for(i in (j+1):p){
    #j = 1, i=2, #3rd column should be to 1
    newX[,i+1] <- ((synX - midpoint)^(i-j))/factorial(i-j)
  }
  #Append newX to X and return
  newX <- rbind(X, newX)
  return(newX)
}

###This buildFeature adds a constant column for me. I do not want glmnet to do it for me
buildFeatureV2 <- function(midpoint, p, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = p+1)
  X[,1] <- rep(1,length(xValues))
  for(i in 1:p){
    X[,i+1] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}


#Testing Purposes
#k <- 18
#deriv1 <- createSyntheticData(sampleData$x, sampleData$y, k)
#plot for testing purposes
#plot(deriv1, col='red')
#lines(gridPoints, trueFunctions[[2]])
#currentPoint <- gridPoints[200]
#X <- buildFeatureV2(currentPoint, p, sampleData$x)

#newX <- addSynthticData(X, deriv1$x, 1, currentPoint)
#newY <- c(sampleData$y, deriv1$y)



k <- 18
synData <- createSyntheticData(sampleData$x, sampleData$y, k)
#synData2 <- createSyntheticDataV2(sampleData$x, sampleData$y, k)
#par(mfrow=c(1,2))
plot(synData, ylim=c(-100,100))
lines(gridPoints, trueFunctions[[2]], lwd=2)
#plot(synData2, ylim=c(-100,100))
#lines(gridPoints, trueFunctions[[2]], lwd=2)


lassoOutput <- matrix(nrow=numGridPoints, ncol=p+1) #ith column shows estimates of (i-1)th derivative 
lambdas <- numeric(numGridPoints) #Tracks lambda parameter over x


#Estimate lasso polynomial at every gridpoint
for(i in 1:numGridPoints){
  
  currentPoint <- gridPoints[i]
  X <- buildFeatureV2(currentPoint, p, sampleData$x)
  X <- addSynthticData(X, synData$x, 1, currentPoint)
  y <- c(sampleData$y, synData$y)
  lassoWeights1 <- computeWeights(sampleData$x, currentPoint, bandwidth, kernel='epak')
  lassoWeights2 <- computeWeights(synData$x, currentPoint, bandwidth, kernel='epak')
  lassoWeights <- c(lassoWeights1, lassoWeights2)
  
  #Penalty factor ensures we dont shrink constant column
  pen_factor <- c(0, rep(1,p))
  #Performing lasso regression, finding optimal lambda
  lassoFit <- glmnet::cv.glmnet(X, y, weights = lassoWeights, maxit=10**7, nlambda=100,
                                intercept=FALSE, penalty.factor=pen_factor)
  lambdas[i] <- lassoFit$lambda.min
  
  # Selects the desired coefficients with the best lambda
  lassoCoef <- coef(lassoFit, s="lambda.min")
  #The i'th row contains all estimated p+1 derivatives estimation at the i'th gridpoint
  lassoOutput[i, ] <- as.vector(lassoCoef)[-1]
}

smoothLambdas <- stats::lowess(gridPoints, lambdas, f=1/10)$y
smoothLambdas <- pmax(rep(0, length(smoothLambdas)), smoothLambdas)
smoothLassoOutput <- matrix(nrow=401, ncol=p+1)

#Refit regression with smoothed lambdas
for(i in 1:length(gridPoints)){
  currentPoint <- gridPoints[i]
  
  currentPoint <- gridPoints[i]
  X <- buildFeatureV2(currentPoint, p, sampleData$x)
  X <- addSynthticData(X, synData$x, 1, currentPoint)
  y <- c(sampleData$y, synData$y)
  lassoWeights1 <- computeWeights(sampleData$x, currentPoint, bandwidth, kernel='epak')
  lassoWeights2 <- computeWeights(synData$x, currentPoint, bandwidth, kernel='epak')
  lassoWeights <- c(lassoWeights1, lassoWeights2)
  
  
  lassoFit <- glmnet::glmnet(X, y, weights = lassoWeights, maxit=10**7, intercept=FALSE, penalty.factor=pen_factor)
  
  # Selects the desired coefficients with the best lambda
  smoothLassoCoef <- coef(lassoFit, s=smoothLambdas[i], exact=TRUE, x=X, y=y, weights=lassoWeights, penalty.factor=pen_factor)
  
  #The i'th row contains all estimated p+1 derivatives estimation at the i'th gridpoint
  smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef)[-1]
}

testComp <- lassoLPR(sampleData$x, sampleData$y, bandwidth)


plotnum <- 1
plot(gridPoints, testComp$lasso[,plotnum], type='l', lwd=2, col='red')
#lines(gridPoints, lassoOutput[,plotnum], type='l', lwd=2, col='blue')
lines(gridPoints, smoothLassoOutput[,plotnum], type='l', lwd=2, col='green')
lines(gridPoints, trueFunctions[[plotnum]], lwd=2)




plot(gridPoints, lambdas, type='l', col='blue')
lines(gridPoints, testComp$lambdas, col='red')
lines(gridPoints, smoothLambdas, col='green')
