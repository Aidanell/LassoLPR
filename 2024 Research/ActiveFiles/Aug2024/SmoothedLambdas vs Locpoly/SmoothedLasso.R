
left <- 0
right <- 1

numTerms <- 10
numPoints = 1000
set.seed(80)
bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2)) #sigma = 0.1
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2)) # sigma = sqrt(0.5)
sine <- expression(sin(5*pi*x)) #sigma = 0.5

equation <- peak
sigma <- sqrt(0.5)

#Calculates all the derivatives needed for plotting exact functions and calculating errors
derivList <- derivCalc(func = equation, numDeriv = numTerms)

##############Building Simulated Data##################
x <- sort(runif(numPoints, min=left, max=right))
noise <- rnorm(n = length(x), mean = 0, sd = sigma)
exactY <- eval(equation)
y <- exactY + noise

locpolyBandwidth <- thumbBw(x,y,0,gaussK)
#lassoBandwidth <- thumbBw(x,y,0,EpaK)
lassoBandwidth <- dpill(x,y)*8
#lassoBandwidth <- locpolyBandwidth
evalPoints <- seq(left, right, length.out=401)

#######Cv.Glmnet calculations########
lassoOutput <- matrix(nrow=401, ncol=numTerms+1)
ListOfLambdas <- c()

for(i in 1:length(evalPoints)){
  currentPoint <- evalPoints[i]
  
  X <- buildFeature(currentPoint, numTerms, x)
  lassoWeights <- computeWeights(x, currentPoint, lassoBandwidth)
  lassoFit <- cv.glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)
  ListOfLambdas <- c(ListOfLambdas, lassoFit$lambda.min)
  
  # Selects the desired coefficients with penalty lambda
  lassoCoef <- coef(lassoFit, s="lambda.min")
  #The i'th row contains all derivative estimation at that specificed evalpoint
  lassoOutput[i, ] <- as.vector(lassoCoef) 
}


##########Computing Smoothed Lambdas Lasso############
smoothedLambas <- lowess(evalPoints, ListOfLambdas, f=1/10)$y
smoothLassoOutput <- matrix(nrow=401, ncol=numTerms+1)

for(i in 1:length(evalPoints)){
  currentPoint <- evalPoints[i]
  
  X <- buildFeature(currentPoint, numTerms, x)
  lassoWeights <- computeWeights(x, currentPoint, lassoBandwidth)
  lassoFit <- glmnet(X, y, weights = lassoWeights, standardize=TRUE, alpha=1, maxit=10**7)

  # Selects the desired coefficients with penalty lambda
  smoothLassoCoef <- coef(lassoFit, s=smoothedLambas[i])
  
  #The i'th row contains all derivative estimation at that specificed evalpoint
  smoothLassoOutput[i, ] <- as.vector(smoothLassoCoef) 
}

plot(evalPoints, lassoOutput[,1])
lines(evalPoints, smoothLassoOutput[,1])
plot(evalPoints, ListOfLambdas)
lines(evalPoints, smoothedLambas)
