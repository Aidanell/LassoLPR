library(locpol)
library(lassoLPR)
buildData <- function(n, equation, sigma){
  x <- sort(runif(n, min=0, max=1))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  trueY <- eval(equation)
  y <- trueY + noise
  return(list(x=x,y=y,trueY=trueY))
}

#Create neccesary data----
n <- 100
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5*x)
p <- 9
sigma <- sqrt(0.5)
deriv <- 3 #Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,EpaK)
gridPoints <- seq(0,1,length.out=401)
#True Functions----
derivCalc <- function(func, numDeriv){
  drvList <- c(func)
  nextDrv <- func
  for(i in 1:numDeriv){
    nextDrv <- D(nextDrv, 'x')
    drvList <- c(drvList, nextDrv)
  }
  return(drvList)
}
derivList <- derivCalc(func = peak, numDeriv = 5)
x <- seq(0,1,length.out=401)
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),eval(derivList[6]))

#Testing----
testLasso <- lassoLPR(sampleData$x, sampleData$y, bandwidth,p)

#Get all local max/min indicies
critIndicies <- which(abs(diff(sign(diff(testLasso$lasso[,4]))))==2)+1
critPoints <- gridPoints[critIndicies]

plot(gridPoints, testLasso$lasso[,4], type='l')
points(critPoints, testLasso$lasso[critIndicies,4], col='red')

#At this point, many of the critical points will be at small jumps that upon visual inspection
#would make it obvious they are not suitable knots. We must filter out small ones via some selection method
#Possible options include top x%, via some threshhold T, or top n points. All sorted based on magnitude of 3rd
#derivative at the critical point


#let's try just the top 20% of crit points and go from there
critValues <- testLasso$lasso[critIndicies,4]
threshhold <- quantile(abs(critValues), 0.4)
knotIndicies <- which(abs(critValues) > threshhold)
knots <- critPoints[knotIndicies]

plot(gridPoints, testLasso$lasso[,4], type='l')
points(knots, critValues[knotIndicies], col='red')

#I also want to get rid of knots that are extremely close to eachother. 
#My current metric will be to take the average of any knots within %1 the range of our data between them.
#This is arbitrary and absolutely should be played around with

range <- diff(range(sampleData$x))
tooClose <- range * 0.01

removeRepeatKnots <- function(knots, tooClose){
  outputKnots <- c()
  while(length(knots) > 0 ){
    currentKnot <- knots[1]
    
    knotDistance <- abs(knots - currentKnot)
    tooCloseKnots <- which(knotDistance < tooClose)
    
    outputKnots <- c(outputKnots, mean(knots[tooCloseKnots])) #Take average of knots within range of eachother 
    
    knots <- knots[! knots %in% knots[tooCloseKnots]] #Removing accounted for knots
  }
  return(outputKnots)
}

finalKnots <- removeRepeatKnots(knots, tooClose)


plot(gridPoints, testLasso$lasso[,4], type='l')
points(finalKnots, rep(0,length(finalKnots)), col='red')

#Now we have our knots needed to fit the natural cubic spline
dataFrame <- data.frame(sampleData)
splineOutput <- lm(y ~ splines::ns(x, knots=finalKnots), data=dataFrame)
plot(sampleData)
lines(gridPoints, predict(splineOutput, data.frame(gridPoints)), col='purple', lwd=2)
points(finalKnots, rep(0,length(finalKnots)), col='red')
