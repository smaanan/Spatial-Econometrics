library(spatstat)
source('ex-5.20-suppl.R')
data(swedishpines)
X <- swedishpines


# Strauss process   
nstep<-70
maxlogpl<-numeric(nstep)
r<-seq(1,15,length=nstep)
for (i in 1:nstep) {
  v<-ppm(X, ~1, interaction=Strauss(r=r[i]),correction='none')
  maxlogpl[i]<-v$maxlogpl     
}   
  rmax<-r[which.max(maxlogpl)]
  rmax<-round(rmax,1)
 plot(r,maxlogpl,type='l', ylab=expression(paste("log PV", (hat(theta)))),xlab="r")
  abline(v=rmax)
    
  fit.mpl<-ppm(X, interaction=Strauss(r=rmax),correction='none') 


  set.seed(30)	
  fit.mcmcmle <- mcmcmle(X, interaction=Strauss(r=rmax),correction='none',nsim=1000,nrmh=1e4,maxiter=1)

  print(fit.mcmcmle)
  fit.ml<-fit.mpl
  fit.ml$fitin$coefs<-fit.mcmcmle$theta
  fit.ml$coef<-fit.mcmcmle$theta
  set.seed(19)	
  env.mcmcmle<-envelope(fit.ml,nsim=40)
  plot(env.mcmcmle,sqrt(cbind(obs,lo,hi)/pi)-r~r,lty=c(1,2,2),ylab=expression(hat(L)(h)-h),xlab='h',col=c(1,1,1), main = ' ')





