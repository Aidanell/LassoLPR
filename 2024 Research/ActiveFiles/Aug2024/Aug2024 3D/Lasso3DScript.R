library("rgl")


#Script Parameters----
gridLength <- 40 #How many grid points in each direction (i.e we will evaluate at gridLength**2 points)
xylim <- c(0,1)
n <- 600
sigma <- sqrt(0.5)
p <- 5 #How many degrees in Lasso

#True Regression Function
f <- function(x,y){return(2-5*x-5*y +5*exp(-20*(x-0.5)**2 -20*(y-0.5)**2))}


#Getting gridPoints and true function data
trueResults <- getTrueFunctionData(xylim, gridLength)
evalPoints <- seq(xylim[1], xylim[2], length.out = gridLength)

evalPointsX <- rep(1:gridLength, gridLength)/gridLength
evalPointsY <- c()
for(i in 1:gridLength){
  evalPointsY <- c(evalPointsY, rep(i,gridLength)/gridLength)
}


#Building data
x <- sort(runif(n, min=xylim[1], max=xylim[2]))
y <- runif(n, min=xylim[1], max=xylim[2])
y <- y[order(match(y, x))]
z <- f(x, y) + rnorm(n, sd=sigma)

generatedData <- data.frame("x"=x, "y"=y, "z"=z)

#Bandwidths
Xbandwidth <- dpill(x,z) * 6
Ybandwidth <- dpill(y,z) * 6

#Building nonparametric Models
lassoResult <- buildLassoMatricies(x,y,z,n,p, Xbandwidth, Ybandwidth, evalPointsX, evalPointsY, gridLength)
unsmoothedLambdas <- lassoResult[[2]]
lassoResult <- lassoResult[[1]]

smoothedLambdas <- loess(unsmoothedLambdas ~ evalPointsX + evalPointsY, span=0.9)$y
smoothLassoResult <- smoothLassoComputation(x,y,z,n,p, smoothedLambdas, Xbandwidth, Ybandwidth, evalPointsX, evalPointsY, gridLength)

#plotting
col <- cm.colors(20)[1 + round(19*(z - min(z))/diff(range(z)))]

open3d()
plot3d(evalPointsX, evalPointsY, as.vector(lassoResult[[1]]), type='s',
       size=1, main="hey", lit=TRUE,sub="3-D Plot")

open3d()
dxyz2 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=lassoResult[[2]])
persp3d(dxyz2, col = col, smooth = FALSE)


open3d()
dxyz3 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=smoothLassoResult[[2]])
persp3d(dxyz3, col = col, smooth = FALSE)


open3d()
dxyz2 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=unsmoothedLambdas)
persp3d(dxyz2, col = col, smooth = FALSE)


open3d()
dxyz3 <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=smoothedLambdas)
persp3d(dxyz3, col = col, smooth = FALSE)



