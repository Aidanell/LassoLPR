
Lasso2DRegression <- function(equation,
                              xylim,
                              numPoints,
                              sigma,
                              bandwidthMulti = 4,
                              seed = as.integer(runif(1,0,1000))
                              ){
  
  set.seed(seed)
  
  exactX <- exactY <- seq(xylim[1], xylim[2], length.out = length)
  exactZ <- outer(exactX,exactY,f)
  evalPoints <- exactX
  
  #Building Data
  x <- sort(runif(numPoints, min=xylim[1], max=xylim[2]))
  y <- runif(numPoints, min=xylim[1], max=xylim[2])
  y <- y[order(match(y, x))]
  z <- f(x, y) + rnorm(numPoints, sd=sigma)
  DataList <- vector(mode='list', length=9)
  
  #Bandwidths
  Xbandwidth <- dpill(x,z) * bandwidthMulti
  Ybandwidth <- dpill(y,z) * bandwidthMulti
  
  #Building Nonparametric models
  lassoModel <- buildLassoMatricies(x, y, z, evalPoints, numPoints, length, Xbandwidth, Ybandwidth)
  domainSize <- xylim[2] - xylim[1]
  NWResult <- buildLocPolyMatricies(x,y,z, degree=3, Xbandwidth/domainSize, Ybandwidth/domainSize)
  
  #Plotting
  par(mfrow=c(2,2))
  persp(exactX, exactY, exactZ, col='lightblue',theta=60, phi=20,
        ticktype='detailed', shade=0.3, main="Exact Graph") -> res
  points(trans3d(x, y, z, pmat = res), col = 2, pch = 16)
  
  persp(evalPoints, evalPoints, lassoModel[[1]], col='lightblue',theta=60, phi=20,
        ticktype='detailed', shade=0.3, main="Lasso Result")
  persp(evalPoints, evalPoints, NWResult[[1]], col='lightblue',theta=60, phi=20,
        ticktype='detailed', shade=0.3, main="Nadaraya-Watson Result")
  
  LassoMSE <- sum( (exactZ - lassoModel[[1]])**2 )
  LocPolyMSE <- sum( (exactZ - NWResult[[1]])**2 )
  
  return(c(LassoMSE, LocPolyMSE))
  
}


f <- function(x,y){return(2-5*x-5*y +5*exp(-20*(x-0.5)**2 -20*(y-0.5)**2))}
length <- 40
xylim <- c(0,1)
numPoints <- 250
sigma <- 0.5

Lasso2DRegression(f, xylim, numPoints, sigma)


repeatSimulations <- function(loops){
  
  ErrorData <- matrix(nrow=loops, ncol = 2)
  for(i in 1:loops){
    nextRun <- Lasso2DRegression(f,xylim, numPoints, sigma)
    print(nextRun)
    ErrorData[i,] <- nextRun
  }
  return(ErrorData)
}

testScript <- repeatSimulations(2)

testScript



bandwidthSim <- function(seed = as.integer(runif(1,0,1000))){
  multipliers <- seq(2,10,1)
  ErrorData <- matrix(nrow=length(multipliers), ncol = 2)
  for(i in 1:length(multipliers)){
    nextRun <- Lasso2DRegression(f,xylim, numPoints, sigma, multipliers[i])
    ErrorData[i,] <- nextRun
  }
  return(ErrorData)
}

#testbandSim <- bandwidthSim()




#plot(seq(2,10,1), testbandSim[,1])
#dev.off()



