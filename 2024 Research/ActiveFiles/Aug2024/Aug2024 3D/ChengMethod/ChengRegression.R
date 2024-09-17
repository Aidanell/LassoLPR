#Testing Cheng's Method


#Script Parameters----
gridLength <- 40 #How many grid points in each direction (i.e we will evaluate at gridLength**2 points)
xylim <- c(0,1)
n <- 250
sigma <- sqrt(0.5)
p <- 5 #How many degrees in Lasso

#True Regression Function
f <- function(x,y){return(2-5*x-5*y +5*exp(-20*(x-0.5)**2 -20*(y-0.5)**2))}

#Getting gridPoints and true function data
evalPointsX <- rep(1:gridLength, gridLength)/gridLength
evalPointsY <- c()
for(i in 1:gridLength){
  evalPointsY <- c(evalPointsY, rep(i,gridLength)/gridLength)
}

trueResults <- getTrueFunctionData(xylim, gridLength, evalPointsX, evalPointsY)
evalPoints <- seq(xylim[1], xylim[2], length.out = gridLength)


#Building data
x <- sort(runif(n, min=xylim[1], max=xylim[2]))
y <- runif(n, min=xylim[1], max=xylim[2])
y <- y[order(match(y, x))]
z <- f(x, y) + rnorm(n, sd=sigma)

generatedData <- data.frame("x"=x, "y"=y, "z"=z)


#Bandwidths
Xbandwidth <- dpill(x,z) * 4
Ybandwidth <- dpill(y,z) * 4

Xbandwidths <- seq(0.1, 1, length.out=20)
Ybandwidths <- seq(0.1, 1, length.out=20)



#Evaluating in the middle of the function
currentX <- evalPointsX[15]
currentY <- evalPointsY[15]



#Build feature matrix and weights
feature <- build2DFeature(currentX, currentY, x, y, p)
#print(feature)

xbands <- c()
ybands <-c()
ghat <- c()
lambdas <- c()

for(i in 1:length(Xbandwidths)){
  for(j in 1:length(Ybandwidths)){
    
    currentXband <- Xbandwidths[i]
    currentYband <- Ybandwidths[j]
    xbands <- c(xbands, currentXband)
    ybands <- c(ybands, currentYband)
    
    
    xWeights <- computeWeights(x, currentX, currentXband)
    yWeights <- computeWeights(y, currentY, currentYband)
    weights <- xWeights * yWeights
    
    #Fit lasso and find best lambda
    set.seed(42)
    lassoFit <- cv.glmnet(feature, z, weights = weights, standardize=TRUE, maxit=1000000)
    lassoCoef <- as.vector(coef(lassoFit, s="lambda.min"))
    
    ghat <- c(ghat, lassoCoef[1])
    lambdas <- c(lambdas, lassoFit$lambda.min)
  }
}

data <- data.frame(x=xbands**(p+1), y=ybands**(p+1), ghat=ghat, lambdas=lambdas)







chengsModel <- lm(ghat ~ x + y, data=data)



library("plot3D")

grid.lines = 40
x.pred <- seq(min(data$x), max(data$x), length.out = grid.lines)
y.pred <- seq(min(data$y), max(data$y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(chengsModel, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)

# create the fitted points for droplines to the surface
fitpoints <- predict(chengsModel)

# scatter plot with regression plane
scatter3D(data$x, data$y, data$ghat, pch = 19, cex = 1,colvar = NULL, col="red", 
          theta =20, phi = 30, bty="b",
          xlab = "X bandwidth", ylab = "Y bandwidth", zlab = "Ghat Estimate",
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = TRUE, fit = fitpoints, col=ramp.col (col = c("dodgerblue3","seagreen2"), n = 300, alpha=0.9), border="black"), main = "Advertising")






