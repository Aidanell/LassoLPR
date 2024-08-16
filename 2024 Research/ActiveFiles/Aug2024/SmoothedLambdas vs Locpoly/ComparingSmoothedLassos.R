left <- 0
right <- 1

numTerms <- 10
numPoints = 100
set.seed(100)
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





#Test 1
#Bandwidth = Dpill * 8, Epan Weights
bandwidth1 <- dpill(x,y)*8
lassoTest1 <- LassoSmoothed(x,y, bandwidth1, TRUE, numTerms)


#Test 2
#Bamdwidth = Thumbbw(degree=0, epaK), Epan Weights
bandwidth2 <- thumbBw(x,y,0,EpaK)
lassoTest2 <- LassoSmoothed(x,y,bandwidth2, TRUE, numTerms)

#Test 3
#Bandwidth = Thumbbw(degree=0, gaussK), Normal Weights
bandwidth3 <- thumbBw(x,y,0,gaussK)
lassoTest3 <- LassoSmoothed(x,y,bandwidth3, FALSE, numTerms)

#Test 4
#Bamdwidth = Thumbbw(degree=3, epaK), Epan Weights
bandwidth4 <- thumbBw(x,y,3,EpaK)
lassoTest4 <- LassoSmoothed(x,y,bandwidth4, TRUE, numTerms)


###Plotting
temp <- x
x <- evalPoints
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]))

x <- temp

deriv <- 0
plot(evalPoints, trueFunctions[[deriv+1]], type='l', lwd=2)
lines(evalPoints, lassoTest1[[2]][,deriv+1], col='blue', lwd=2)
lines(evalPoints, lassoTest2[[2]][,deriv+1], col='orange', lwd=2)
#lines(evalPoints, lassoTest2[[2]][,deriv+1], col='red', lwd=2)
lines(evalPoints, lassoTest3[[2]][,deriv+1], col='red', lwd=2)

