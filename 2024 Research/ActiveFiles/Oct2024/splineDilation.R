#Create neccesary data----
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5*x)
p <- 9
sigma <- sqrt(0.5)
deriv <- 3 #Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,EpaK)
gridPoints <- seq(0,1,length.out=401)

testLasso <- lassoLPR(sampleData$x, sampleData$y, bandwidth,p)

#Get all local max/min indicies
critIndicies <- which(abs(diff(sign(diff(testLasso$lasso[,4]))))==2)+1
critPoints <- gridPoints[critIndicies]

plot(gridPoints, testLasso$lasso[,4], type='l')
points(critPoints, testLasso$lasso[critIndicies,4], col='red')

maxDilation <- VPdtw::dilation(testLasso$lasso[,4], 5)
minDilation <- VPdtw::dilation(-testLasso$lasso[,4], 5)
dilationInvariant <- which(testLasso$lasso[,4] == maxDilation | -testLasso$lasso[,4] == minDilation)

#Only points that are critical and dilation invariant
knotIndex <- critIndicies[which(critIndicies %in% dilationInvariant)]
knots <- gridPoints[knotIndex]
knotSize <- testLasso$lasso[knotIndex,4]

plot(gridPoints, testLasso$lasso[,4], type='l')
points(knots, testLasso$lasso[knotIndex,4], col='red')

threshhold <- quantile(abs(critValues), 0.5)
knotIndex <- knotIndex[which(abs(knotSize) > threshhold)]
knots <- gridPoints[knotIndex]


dataFrame <- data.frame(sampleData)
splineOutput <- lm(y ~ splines::ns(x, knots=knots), data=dataFrame)
splinePred <- predict(splineOutput, data.frame(gridPoints))
plot(sampleData)
lines(gridPoints, splinePred, col='purple', lwd=2)
points(knots, splinePred[knotIndex], col='red')



