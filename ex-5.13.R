library(spdep)
dyn.load(paste("ex-5.13-suppl",.Platform$dynlib.ext,sep=""))
rm(list=ls())
mercer.hall<-read.table('mercer-hall.txt',header=FALSE)
names(mercer.hall)<-c("x","y","z")
n<-20
m<-25
z<-mercer.hall$z
Z<-t(matrix(z,m,n))
zd<-Z/max(Z)
len<-2.0*(mercer.hall$z-min(mercer.hall$z))/(max(mercer.hall$z)-min(mercer.hall$z))
plot(mercer.hall$x,mercer.hall$y,type='n',xlab="colonne",ylab="ligne",axes=FALSE,frame.plot=TRUE)
symbols(mercer.hall$x,mercer.hall$y,squares=len,xlab="colonne",ylab="ligne",inches=FALSE, add = TRUE)


plot(mercer.hall$x,mercer.hall$y,type='n',xlab="colonne",ylab="ligne",axes=FALSE,frame.plot=TRUE)
coding<-seq(1,500,2)
points(mercer.hall$x[coding],mercer.hall$y[coding],pch=3)
coding<-seq(2,500,2)
points(mercer.hall$x[coding],mercer.hall$y[coding],pch=20)


plot(z~as.factor(x),data=mercer.hall,axes=FALSE,frame.plot=TRUE,xlab='colonne')
axis(1, at=1:m,labels = 1:m)


plot(z~as.factor(y),data=mercer.hall,axes=FALSE,frame.plot=TRUE,xlab="ligne",horizontal=TRUE)
axis(1, at=1:n,labels = 1:n)
fit.lm<-lm(z~as.factor(x),data=mercer.hall)

z<-fit.lm$resid
Z<-t(matrix(z,m,n))

neighborhood.type<-0
nparam<-2
edge<-1
border<-2*edge


x<-numeric(n*m*nparam)
# toroidal observations
ymat<-matrix(0,n+border,m+border)
ymat[2:(n+1),2:(m+1)]<-Z
ymat[1,2:(m+1)]<-Z[n,]
ymat[n+border,2:(m+1)]<-Z[1,]
ymat[2:(n+1),1]<-Z[,m]
ymat[2:(n+1),m+border]<-Z[,1]
mu.hat<-mean(Z)
c00<-mean(Z^2)-(mu.hat)^2


n<-n+border
m<-m+border

xdata<-.C("matneigh", ymat=as.double(ymat),x=as.double(x),n=as.integer(n), m=as.integer(m),neighborhood.type=as.integer(neighborhood.type))
	n<-n-border
	m<-m-border
	x<-matrix(xdata$x,n*m,nparam)
y<-as.vector(t(ymat[(edge+1):(n+edge),(edge+1):(m+edge)]))


s.eig<-function(omega,eta,theta)
{
a<-2*theta[1]*(cos(omega))+2*theta[2]*(cos(eta))
return(a)
}    

C.eig<-function(n,m,theta) {
    a<-numeric(n*m)
    omega<-2*pi*(0:(n-1))/n
    eta<-2*pi*(0:(m-1))/m 
    a<-1-outer(omega,eta,s.eig,theta)
    return(a)       
}   

mu.hat<-mean(Z)

c10<-0
c01<-0
for (i in 1:n) 
   for (j in 1:m)
   {

    c10<-c10+(ymat[i+edge,j+edge]-mu.hat)*(ymat[(i+1+edge),(j+edge)]-mu.hat)
    c01<-c01+(ymat[i+edge,j+edge]-mu.hat)*(ymat[(i+edge),(j+edge+1)]-mu.hat)
   }
c10<-c10/(n*m)
c01<-c01/(n*m)   
nobs<-n*m
neg.loglik<-function(theta){


      b<-C.eig(n,m,theta)
      if( sum(b < 0) >0 ) 
      {
      
      }
      b<-prod(b)
      temp<-c00-2*(theta[1]*c10+theta[2]*c01)

      if( temp <=0 ) 
      {
      return(10e300)
      }
      val<-log(temp)-log(b)/nobs
      return(val)
}

fit.mpl.0<-lm(y~-1)
fit.cod1.0<-lm(y~-1, subset=seq(1,500,2))
fit.cod2.0<-lm(y~-1, subset=seq(2,500,2))
u.mpl.0<-logLik(fit.mpl.0)
u.cod1.0<-logLik(fit.cod1.0)
u.cod2.0<-logLik(fit.cod2.0)
u.mle.0<-(-0.5*n*m)*log(c00)
fit.mpl<-lm(y~x-1)
fit.cod1<-lm(y~x-1, subset=seq(1,500,2))
fit.cod2<-lm(y~x-1, subset=seq(2,500,2))
u.mpl<-logLik(fit.mpl)
u.cod1<-logLik(fit.cod1)
u.cod2<-logLik(fit.cod2)
sigma2.mpl<-sqrt(mean(fit.mpl$resid^2))
sigma2.cod1<-sqrt(mean(fit.cod1$resid^2))
sigma2.cod2<-sqrt(mean(fit.cod2$resid^2))
theta.mpl<-fit.mpl$coef[1:2]
theta.cod1<-fit.cod1$coef[1:2]
theta.cod2<-fit.cod2$coef[1:2]
theta<-c(fit.mpl$coef[1:2])
a<-optim(theta,neg.loglik,method =  "BFGS")
theta.mle<-a$par
sigma2.mle<-c00-2*(theta.mle[1]*c10+theta.mle[2]*c01)
u.mle<-(-0.5*n*m)*a$value


x<-x[,1]+x[,2]
fit.mpl.iso<-lm(y~x-1)
fit.cod1.iso<-lm(y~x-1, subset=seq(1,500,2))
fit.cod2.iso<-lm(y~x-1, subset=seq(2,500,2))
theta.mpl.iso<-fit.mpl.iso$coef[1]
theta.cod1.iso<-fit.cod1.iso$coef[1]
theta.cod2.iso<-fit.cod2.iso$coef[1]

u.mpl.iso<- as.numeric(logLik(fit.mpl.iso))
u.cod1.iso<- as.numeric(logLik(fit.cod1.iso))
u.cod2.iso<- as.numeric(logLik(fit.cod2.iso))

s.eig<-function(omega,eta,theta)
{
a<-2*theta[1]*(cos(omega)+cos(eta))
return(a)
}    

neg.loglik<-function(theta){


      b<-C.eig(n,m,theta)
      if( sum(b < 0) >0 ) 
      {

      }
      b<-prod(b)
      temp<-c00-2*(theta[1]*c10+theta[1]*c01)

      if( temp <=0 ) 
      {

      return(10e300)
      }
      val<-log(temp)-log(b)/nobs
      return(val)
}
theta<-c(fit.mpl.iso$coef[1:1])
a<-optim(theta,neg.loglik,method =  "BFGS")
theta.mle.iso<-a$par
sigma2.mle.iso<-c00-2*(theta.mle[1]*c10+theta.mle[1]*c01)
u.mle.iso<-(-0.5*n*m)*a$value

results<-c(NA,NA,u.cod1.0,NA,NA,u.mle.0)
results<-rbind(results,c(theta.cod1.iso,NA,u.cod1.iso,theta.mle.iso,NA,u.mle.iso))
results<-rbind(results,c(theta.cod1,u.cod1,theta.mle,u.mle))


results<-as.numeric(c(NA,NA,u.cod1.0,NA,NA,u.cod2.0,NA,NA,u.mle.0))
results<-rbind(results,as.numeric(c(theta.cod1.iso,NA,u.cod1.iso,theta.cod2.iso,NA,u.cod2.iso,theta.mle.iso,NA,u.mle.iso)))
results<-rbind(results,as.numeric(c(theta.cod1,u.cod1,theta.cod2,u.cod2,theta.mle,u.mle)))


