

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


t <- 40
possibleBandwidths <- seq(bandwidth*0.1, bandwidth*3, length.out=t)
lassoOutput <- matrix(nrow=t, ncol=p+1)
lambdas <- numeric(t)
for(i in 1:t){
  currBand <- possibleBandwidths[i]
  
  currentPoint <- gridPoints[25]
  
  X <- buildFeature(currentPoint, p, sampleData$x)
  lassoWeights <- computeWeights(sampleData$x, currentPoint, currBand, "epak")
  
  #Performing lasso regression, finding optimal lambda
  set.seed(10)
  lassoFit <- glmnet::cv.glmnet(X, sampleData$y, weights = lassoWeights, maxit=10**7, nlambda=100)
  lambdas[i] <- lassoFit$lambda.min
  
  # Selects the desired coefficients with the best lambda
  lassoCoef <- coef(lassoFit, s="lambda.min")
  #The i'th row contains all estimated p+1 derivatives estimation at the i'th gridpoint
  lassoOutput[i, ] <- as.vector(lassoCoef) 
  
}

test <- 2
plot(possibleBandwidths, lassoOutput[,test], type='l')
abline(trueFunctions[[test]][25],0, col='red')
abline(v=bandwidth, col='blue')


plot(possibleBandwidths, lambdas, type='l')
abline(v=bandwidth, col='blue')
