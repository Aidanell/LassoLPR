library("rgl")
library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)
library(locpol)

#Script Parameters----
gridLength <- 40 #How many grid points in each direction (i.e we will evaluate at gridLength**2 points)
xylim <- c(0,1)
n <- 250
sigma <- sqrt(0.5)
p <- 10 #How many degrees in Lasso

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

#Building nonparametric Models
lassoResult <- buildLassoMatricies(x,y,z,n,p, Xbandwidth, Ybandwidth, evalPointsX, evalPointsY, gridLength)
unsmoothedLambdas <- lassoResult[[2]]
lassoResult <- lassoResult[[1]]

smoothedLambdas <- loess(unsmoothedLambdas ~ evalPointsX + evalPointsY, span=0.1)$fitted
smoothLassoResult <- smoothLassoComputation(x,y,z,n,p, smoothedLambdas, Xbandwidth, Ybandwidth, evalPointsX, evalPointsY, gridLength)

#plotting
col <- cm.colors(20)[1 + round(19*(z - min(z))/diff(range(z)))]
deriv <- 5

open3d()
dxyz1 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=lassoResult[[deriv]])
persp3d(dxyz1, col = col, smooth = FALSE, main="Lasso Unsmoothed Lambdas")


open3d()
dxyz2 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=smoothLassoResult[[deriv]])
persp3d(dxyz2, col = col, smooth = FALSE, main="Lasso Smoothed Lambdas")

open3d()
dxyz3 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=trueResults[[deriv]])
persp3d(dxyz3, col = col, smooth = FALSE)


open3d()
dxyz4 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=unsmoothedLambdas)
persp3d(dxyz4, col = col, smooth = FALSE, main="Unsmoothed Lambdas")


open3d()
dxyz5 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=smoothedLambdas)
persp3d(dxyz5, col = col, smooth = FALSE, main="Smoothed Lambdas")



#Locpoly 3D
LLResult <- interp::locpoly(x,y,z, degree=3, pd='all', h=c(Xbandwidth, Ybandwidth))
open3d()
dxyz6 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=as.vector(LLResult$zxy))
persp3d(dxyz6, col = col, smooth = FALSE, main="Locpoly", axes=TRUE, box=TRUE)
