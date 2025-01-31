

##################Helper Methods###########################----
buildFeature <- function(midpoint, p, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = p)
  for(i in 1:p){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}


buildPolynomial <- function(x, coefs, gridPoint){
  output <- coefs[1]
  for(i in 2:length(coefs)){
    output <- output + coefs[i]*(x-gridPoint)^i / factorial(i)
  }
  return(output)
}
#curve(buildPolynomial(x, coefs=c(1,2,3), gridPoint=2), from=0, to=4)


#Function of the Epanechnikov kernel
epan <- function(x, bandwidth){
  ifelse(abs(x) <= 1, 3/4 * (1 - x^2), 0)
}


#Computes the scaled kernel weights for every input x
#Works for Epanechnikov or normal kernels
computeWeights <- function(x, x0, bandwidth, kernel='epak'){
  weights <- numeric(length(x))
  
  for(i in 1:length(x)){
    diff <- abs(x[i]-x0) / bandwidth
    if(kernel == 'epak'){
      weights[i] <- epan(diff)/bandwidth
    }else{
      weights[i] <- dnorm(diff)/bandwidth
    }
  }
  return(weights)
}


#Function that builds simulated data
buildData <- function(n, equation, sigma){
  x <- sort(runif(n, min=0, max=1))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  trueY <- eval(equation)
  y <- trueY + noise
  return(list(x=x,y=y,trueY=trueY))
}


#True Functions
derivCalc <- function(func, numDeriv){
  drvList <- c(func)
  nextDrv <- func
  for(i in 1:numDeriv){
    nextDrv <- D(nextDrv, 'x')
    drvList <- c(drvList, nextDrv)
  }
  return(drvList)
}

############################CV Methods#########################----
AidansCV <-function(x, y, gridPoint, lambdas, weights, p, trueY){
  
  
  #lambdas <- seq(2,0,length.out=20)
  AVEWCV <- numeric(20)
  X <- buildFeature(gridPoint, p, x) #X does not change, we just have to leave out one row at a time
  trueY <- unlist(trueY)
  
  
  
  for(lambda in lambdas){
    WMSE <- numeric(length(y))
    for(i in 1:length(y)){
      
      currentLeftOut <- y[i]
      
      lassoFit <- glmnet::glmnet(X[-i,], y[-i], weights = weights[-1], maxit=10**7)
      lassoCoef <- as.vector(coef(lassoFit, s=lambda))
      evaluateModel <- buildPolynomial(x[-1], lassoCoef, gridPoint)
      
      #print(length(trueY[-1]))
      #print(length(weights[-1]))
      #print(length(evaluateModel))
      
      WMSE[i] <- sum(weights[-1]*(trueY[-1] - evaluateModel)**2)
    }
    #print(lambda)
    AVEWCV <- mean(WMSE)
  }
  
  #Specify Best Lambda
  lambdaIndex <- which.min(AVEWCV)
  chosenLambda <- lambdas[lambdaIndex]
  print(chosenLambda)
  #Use lambda on full dataset
  lassoFit <- glmnet::glmnet(X, y, weights = weights, maxit=10**7)
  lassoCoef <- coef(lassoFit, s=lambda)
  
  
  return(list(lassoCoef, chosenLambda))
}


testingCustomCV <- function(x, y, h, p=10, numGridPoints=401, varBands = FALSE){
  
  gridPoints <- seq(min(x), max(x), length.out=numGridPoints)
  
  lassoOutput <- matrix(nrow=numGridPoints, ncol=p+1) #ith column shows estimates of (i-1)th derivative 
  lambdas <- numeric(numGridPoints) #Tracks lambda parameter over x
  
  
  #Estimate lasso polynomial at every gridpoint
  for(i in 1:numGridPoints){
    
    currentPoint <- gridPoints[i]
    X <- buildFeature(currentPoint, p, x)
    if(varBands == TRUE){
      variableH <- h*variableBands[i]
      lassoWeights <- computeWeights(x, currentPoint, variableH, kernel='norm')
    }else{
      lassoWeights <- computeWeights(x, currentPoint, h, kernel='norm') 
    }
    
    #It does not matter how well the points do if they are well outside the neighborhood. Still experimenting with this
    closeIndex <- which(lassoWeights > 0)
    #closeIndex <- which(lassoWeights > quantile(lassoWeights, 0.5)) #I do not care about performance if they are not in the window
    #closeIndex <- which( abs(x - currentPoint)/1.3*h < 1)
    
    #Calculate the sequence of lambdas to use
    lambdaSeq <- glmnet::glmnet(X[closeIndex,], y[closeIndex], weights = lassoWeights[closeIndex], maxit=10**7, nlambda=25)$lambda
   
    #Perform n-fold cross validation to calculate the 
    PerformCV <- NFoldCV_V2(x[closeIndex], y[closeIndex], currentPoint, lambdaSeq, lassoWeights[closeIndex], p, trueY[closeIndex], folds=5)
    lassoCoef <- PerformCV[[1]]
    lambdas[i] <- PerformCV[[2]]
    print(i)

    
    #The i'th row contains all estimated p+1 derivatives estimation at the i'th gridpoint
    lassoOutput[i, ] <- as.vector(lassoCoef) 
  }
  
  smoothLambdas <- stats::lowess(gridPoints, lambdas, f=1/10)$y
  smoothLassoOutput <- matrix(nrow=401, ncol=p+1)
  
  #Refit regression with smoothed lambdas
  for(i in 1:length(gridPoints)){
    currentPoint <- gridPoints[i]
    
    currentPoint <- gridPoints[i]
    X <- buildFeature(currentPoint, p, x)
    if(varBands == TRUE){
      variableH <- h*variableBands[i]
      lassoWeights <- computeWeights(x, currentPoint, variableH, kernel='norm')
    }else{
      lassoWeights <- computeWeights(x, currentPoint, h, kernel='norm') 
    }
    
    lassoFit <- glmnet::glmnet(X, y, weights = lassoWeights, maxit=10**7)
    #lassoCoef <- coef(lassoFit, s=lambda, exact=TRUE, x=X, y=y, weights=weights)
    
    # Selects the desired coefficients with the best lambda
    smoothLassoCoef <- coef(lassoFit, s=smoothLambdas[i], exact=TRUE, x=X, y=y, weights=lassoWeights)
    
    #The i'th row contains all estimated p+1 derivatives estimation at the i'th gridpoint
    smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef) 
  }
  
  return(list(lassoOutput, lambdas, smoothLassoOutput, smoothLambdas))
}


NfoldCV <-function(x, y, gridPoint, lambdas, weights, p, trueY, folds = 10){
  
  
  #Split data into n folds
  splitData <- split(1:length(x), sample(1:folds, length(x), replace=TRUE))
  
  #Demonstrates how weights are dispoportionate
  #totalWeight <- numeric(folds)
  #for(i in 1:folds){
  #  totalWeight[i] <- sum(weights[splitData[[i]]])
  #}
  #barplot(totalWeight)
  
  AVEWCV <- c()
  X <- buildFeature(gridPoint, p, x) #X does not change, we just have to leave out one row at a time
  
  
  for(lambda in lambdas){
    WMSE <- numeric(folds)
    
    for(i in 1:folds){
      
      currentLeftOut <- as.vector(splitData[[i]])
      
      #Fit model without one of the folds
      lassoFit <- glmnet::glmnet(X[-currentLeftOut,], y[-currentLeftOut], weights = weights[-currentLeftOut], maxit=10**7)
      lassoCoef <- as.vector(coef(lassoFit, s=lambda))
      #Evaluate model with points that were not trained with
      evaluateModel <- buildPolynomial(x[currentLeftOut], lassoCoef, gridPoint)
      
      #Get weighted MSE (with more weight towards points closer to gridPoint)
      WMSE[i] <- sum(weights[currentLeftOut]*(trueY[currentLeftOut] - evaluateModel)**2)
    }
    #Add CVMSE to list
    AVEWCV <- c(AVEWCV, mean(WMSE))
  }
  
  #Specify Best Lambda
  lambdaIndex <- which.min(AVEWCV)
  chosenLambda <- lambdas[lambdaIndex]
  
  print(length(lambdas))
  plot(1:length(lambdas),lambdas)
  points(lambdaIndex, chosenLambda, col='red', cex=3, pch=16)
  
  #Use lambda on full dataset
  lassoFit <- glmnet::glmnet(X, y, weights = weights, maxit=10**7)
  lassoCoef <- coef(lassoFit, s=lambda, exact=TRUE, x=X, y=y, weights=weights)
  
  
  return(list(lassoCoef, chosenLambda))
}






