rm(list=ls())
library(spatstat)
data(betacells)
X <- betacells
plot(X$x,X$y, xlab='x',ylab='y',type='n', xlim=X$window$xrange, ylim= X$window$yrange)
points(X$x[X$marks =='on'],X$y[X$marks =='on'], pch='+',cex=1.2)
points(X$x[X$marks =='off'],X$y[X$marks =='off'], pch='-',cex=1.2)


X <- unmark(betacells)
fit.0<-ppm(X,~1, interaction=Poisson())
fit.1<-ppm(X,~1+x, interaction=Poisson())


vvv<-diagnose.ppm(fit.0,type = "pearson",which='marks',plot.it=FALSE)


par(mar = c(3,3,1,1))
#plot(vvv,srange=range(vvv$Ydens$v),monochrome=TRUE,main=' ')
plot(vvv,main=' ',monochrome=TRUE)


vvv<-diagnose.ppm(fit.0,type = "raw",which='smooth',plot.it=FALSE)


par(mar = c(3,3,1,1))
plot(vvv,main=' ',monochrome=TRUE)


vvv<-diagnose.ppm(fit.1,type = "pearson",which='marks',plot.it=FALSE)


par(mar = c(3,3,1,1))
plot(vvv,monochrome=TRUE,main=' ',srange=range(vvv$Ydens$v,na.rm=TRUE))


vvv<-diagnose.ppm(fit.1,type = "raw",which='smooth',plot.it=FALSE)


par(mar = c(3,3,1,1))
plot(vvv,main=' ',monochrome=TRUE)


lurking(fit.0,expression(x),type = "pearson")


lurking(fit.0,expression(y),type = "pearson")


lurking(fit.1,expression(x),type = "pearson")


lurking(fit.1,expression(y),type = "pearson")


vvv<-qqplot.ppm(fit.1,type = "pearson",verbose=FALSE,plot.it=FALSE)


plot(vvv,pch=20,col=1)


# Strauss process   
nstep<-100
maxlogpl<-numeric(nstep)
r<-seq(1,100,length=nstep)
rmax<-74


fit.mpl<-ppm(X,~1+x, interaction=Strauss(r=rmax),correction = "isotropic") 
vvv<-qqplot.ppm(fit.mpl,type = "pearson",verbose=FALSE,plot.it=FALSE)


plot(vvv,pch=20,col=1)


vvv<-diagnose.ppm(fit.mpl,type = "pearson",which='smooth',plot.it=FALSE)


par(mar = c(3,3,1,1))
plot(vvv,main=' ',monochrome=TRUE)


vvv<-diagnose.ppm(fit.mpl,type = "pearson",which='marks',plot.it=FALSE)


par(mar = c(3,3,1,1))
plot(vvv,monochrome=TRUE,main=' ',srange=range(vvv$Ydens$v,na.rm=TRUE))


