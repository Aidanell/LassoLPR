x <- runif(1000, min=0,1)
x <- rnorm(1000, 0.5, 0.3)
x0 = 0.5
bandwidth = 0.3

gaussKernel <- function(x){
  1/sqrt(2*pi) * exp(-x^2 / 2)
}

Snj <- function(j){
  
  output = sum((x-x0)**j * gaussKernel((x-x0)/bandwidth))
  return(output)
}

buildFeature <- function(midpoint, p, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = p)
  for(i in 1:p){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}

computeWeights <- function(x, x0, bandwidth, kernel='epak'){
  weights <- numeric(length(x))
  
  for(i in 1:length(x)){
    diff <- abs(x[i]-x0) / bandwidth
    if(kernel == 'epak'){
      weights[i] <- epan(diff)/bandwidth
    }else{
      weights[i] <- dnorm(diff)/bandwidth
    }
  }
  return(weights)
}



p <- 9
Sn <- matrix(nrow = p+1, ncol=p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    Sn[i,j] <- Snj(i+j-2)
  }
}

#Remove 3rd and 4th derivatives (as example)
termsToRemove <- c(3,4,7)
SnStar <- Sn[-termsToRemove, -termsToRemove]

#Inverse
SnStarInv <- solve(SnStar)
SnInv <- solve(Sn)

#For Calculating Variance
X <- buildFeature(x0, p, x)
X <- cbind(rep(1,n), X)
XStar <- X[, -termsToRemove]
W <- diag(computeWeights(x, x0, bandwidth, kernel='norm'))

#Asymptotic Variance
VarStar <- SnStarInv %*% t(XStar) %*% W %*% W %*% XStar %*% SnStarInv
Var <- SnInv %*% t(X) %*% W %*% W %*% X %*% SnInv
print(VarStar)
print(Var)


varBetaHats <- diag(Var)
varBetaHatsStar <- diag(VarStar)
#Adds zero where dropped columns are. For easier comparison between the two
for(i in 1:length(termsToRemove)){
  currentTerm <- termsToRemove[i]
  varBetaHatsStar <- append(varBetaHatsStar, 0, after=(currentTerm-1))
}
print(varBetaHatsStar)
print(varBetaHats)

#Now lets calculate bias

largeMoments <- matrix(nrow=(p+1), ncol=1)
for(i in 1:(p+1)){
  largeMoments[i,1] <- Snj(p+i)
}

bias <- as.vector(SnInv %*% largeMoments)

largeMomentsStar <- largeMoments[-termsToRemove, 1]
biasStar <- as.vector(SnStarInv %*% largeMomentsStar)
for(i in 1:length(termsToRemove)){
  currentTerm <- termsToRemove[i]
  biasStar <- append(biasStar, 0, after=(currentTerm-1))
}


