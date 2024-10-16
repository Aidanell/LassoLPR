h
#p=1, r=0
h1 <- thumbBw(x,y,deg=1, gaussK)

n <- length(x)

#Asymptotics say h3 = h1 * (n/(2r+1))^(2/(2p+3)(2p+5)) for p-r odd
#Asymptotics say h3 = h1 8 (n/(2r+1))^(4/(2p+5)(2p+9)) for p-r even

#h3 = h1 * n^(2/35)
#h3 = h1 * n^(2/25)

#h5 = h3* n^(2/99)
#h5 = h3 * n^(4/165)

h3Asymptotics <- h1 * n**(4/117)
h3 <- thumbBw(x,y,deg=3, gaussK)

#h1locpoly <- KernSmooth::locpoly(x,y,degree=1,drv=0, bandwidth=h3, gridsize=401)
h3Asymlocpoly <- KernSmooth::locpoly(x,y,degree=3,drv=0, bandwidth=h3Asymptotics, gridsize=401)
h3locpoly <- KernSmooth::locpoly(x,y,degree=3,drv=0, bandwidth=h3, gridsize=401)



#h5Asymptotics <- h3 * n**(4/165)
h5Asymptotics <- h3 * n**(2/99)
h5 <- thumbBw(x,y,deg=5, gaussK)

  
h5locpolyAsym <- KernSmooth::locpoly(x,y,degree=5,drv=0, bandwidth=h5Asymptotics, gridsize=401)
h5locpoly <- KernSmooth::locpoly(x,y,degree=5,drv=0, bandwidth=h5,gridsize=401)


temp <- x
x <- evalPoints
plot(evalPoints, eval(derivList[deriv+1]), type='l', lwd=2)


#lines(h3Asymlocpoly, type='l', col="purple")
#lines(h3locpoly, type='l', col="red")

lines(h5locpolyAsym, type='l', col="purple")
lines(h5locpoly, type='l', col="red")


x <- temp
