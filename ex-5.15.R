par(mar = c(3,3,1,1))
library(spatstat)
data(bei)
plot(bei$x,bei$y,pch='.',xlab='x',ylab='y')


par(mar = c(3,3,1,1))
  sigma <- (1/8) * c(diff(bei$window$xrange),diff(bei$window$yrange))
  print(round(sigma,2))
  bei.smooth <- density.ppp(bei,sigma)
contour(bei.smooth,main="Noyau anisotropique")
points(bei,pch='.')


library(fields)
elevation<-bei.extra$elev
col<-gray((4:31)/31)
 image.plot(elevation$xcol,elevation$yrow,t(elevation$v),xlab='x',ylab='y',col =col)



gradient<-bei.extra$grad
image.plot(gradient$xcol,gradient$yrow,t(gradient$v),xlab='x',ylab='y',col =col)


bei.fit<-ppm(bei,~elevation+gradient,covariates=list(elevation=elevation,gradient=gradient))
summary(bei.fit)
fisher.inf<-vcov(bei.fit) #Fisher information matrix
sqrt(diag(fisher.inf)) #standard errors
rho<-predict.ppm(bei.fit)
nsim<- 39
env<-envelope(bei, Kinhom,nsim=nsim, lambda=rho, simulate=expression(rpoispp(rho)))

plot(env,cbind(obs,lo,hi)~r,ylab="K(h)",xlab='h', main= ' ', lty=c(1,2,2),col=c(1,1,1))


