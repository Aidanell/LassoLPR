#Script Parameters----
length <- 40
xylim <- c(0,1)
numPoints <- 250
sigma <- sqrt(0.5)

#Defining matricies of true data for the function and all its derivatives----
trueX <- trueY <- seq(xylim[1], xylim[2], length.out = length)
trueResult <- vector(mode='list', length=10)
trueResult[[1]] <- outer(trueX, trueY, f)
trueResult[[2]] <- outer(trueX, trueY, fx)
trueResult[[3]] <- outer(trueX, trueY, fy)
trueResult[[4]] <- outer(trueX, trueY, fxx)
trueResult[[5]] <- outer(trueX, trueY, fxy)
trueResult[[6]] <- outer(trueX, trueY, fyy)
trueResult[[7]] <- outer(trueX, trueY, fxxx)
trueResult[[8]] <- outer(trueX, trueY, fxxy)
trueResult[[9]] <- outer(trueX, trueY, fxyy)
trueResult[[10]] <- outer(trueX, trueY, fyyy)

#Script----
#x/y values we evaluate the function at
evalPoints <- trueX

#Building Data
x <- sort(runif(numPoints, min=xylim[1], max=xylim[2]))
y <- runif(numPoints, min=xylim[1], max=xylim[2])
y <- y[order(match(y, x))]
z <- f(x, y) + rnorm(numPoints, sd=sigma)

generatedData <- data.frame("x"=x, "y"=y, "z"=z)


#Bandwidths
Xbandwidth <- dpill(x,z) * 6
Ybandwidth <- dpill(y,z) * 6



#Building Nonparametric models

lassoResult <- buildLassoMatricies(x, y, z, evalPoints, numPoints, length, Xbandwidth, Ybandwidth)
Lambdas3D <- lassoResult[[2]]
lassoResult <- lassoResult[[1]]
domainSize <- xylim[2] - xylim[1]
locpolyResult <- buildLocPolyMatricies(x,y,z, degree=3, Xbandwidth/domainSize, Ybandwidth/domainSize)
#ErrorDf <- CalculateErrors(trueResult, lassoResult, locpolyResult)

#Creating object that holds all data from simulation
LPRResults <- vector(mode="list", 0)
LPRResults$generatedData <- generatedData
#LPRResults$ErrorDf <- ErrorDf
LPRResults$lassoResult <- lassoResult
LPRResults$locpolyResult <- locpolyResult

#Defining xylim's to better see plots
flim <- c(-8,2.1)
fxlim <- c(-25, 15)
fylim <- c(-25,15)
fxxxlim <- c(-1500,1500)


smoothedLambas <- loess(as.vector(Lambdas3D) ~ evalPointsX + evalPointsY, span=0.1)$y


persp(evalPoints, evalPoints, lassoResult[[1]], col='lightblue',theta=30, phi=20,
      ticktype='detailed', shade=0.3, main="Lasso Result")


persp(evalPoints, evalPoints, Lambdas3D, col='lightblue',theta=30, phi=20,
      ticktype='detailed', shade=0.3, main="Lambdas")


persp(evalPoints, evalPoints, smoothedLambas, col='lightblue',theta=30, phi=20,
      ticktype='detailed', shade=0.3, main="Lambdas")


evalPointsX <- rep(1:40, 40)/40
evalPointsY <- c()
for(i in 1:40){
  evalPointsY <- c(evalPointsY, rep(i,40)/40)
}


library("rgl")

plot3d(evalPointsX, evalPointsY, as.vector(lassoResult[[2]]), type='s',
       size=1, main="hey", lit=TRUE,sub="3-D Plot")


plot3d(evalPointsX, evalPointsY, as.vector(smoothedLambas), type="s",  
       size=1, main="Lambda Magnitude over surface", lit=TRUE,sub="3-D Plot", zlab="Lambda Magnitude")

plot3d(evalPointsX, evalPointsY, as.vector(Lambdas3D), type="s",  
       size=1, main="Lambda Magnitude over surface", lit=TRUE,sub="3-D Plot", zlab="Lambda Magnitude")



plot3d(evalPointsX, evalPointsY, fxx(evalPointsX, evalPointsY), type="s",  
       size=1, main="hey", lit=TRUE,sub="3-D Plot", zlab="Lambda Magnitude")

#s <- scene3d()
#s$par3d$windowRect <- 2*s$par3d$windowRect
#plot3d(s)

open3d()
col <- cm.colors(20)[1 + round(19*(z - min(z))/diff(range(z)))]
dxyz <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=as.vector(lassoResult[[1]]))
persp3d(dxyz, col = col, smooth = FALSE)


dxyz <- deldir::deldir(x=evalPointsX, y=evalPointsY, z=f(evalPointsX, evalPointsY))
persp3d(dxyz, col = col, smooth = FALSE)




