x <- c(0,1,2)
y <- c(4,5,6)

Betas <- lm(y~x)
simga <- 2
Betas$coefficients[1]

Beta1ols <- Betas$coefficients[2]

#MSE(Beta0 Ridge)
curve(sigma*(3*x^2 + 12*x + 30)/(3*x + 6)^2 + (3*x/(6+3*x) * Beta1ols)^2, from=-0.1, to=0.1)
Betas$fitted.values

MSE <- sum(Betas$fitted.values - y)
sigmaXSq <- sum(x**2) - mean(x)^2


Beta1RidgeCurve <- function(x){
                    (MSE*sigmaXSq)/(sigmaXSq + x)^2 +
                    (Beta1ols * (sigmaXSq/(sigmaXSq + x) - 1))^2
}

curve(Beta1RidgeCurve, from=-.01, to=.01)
