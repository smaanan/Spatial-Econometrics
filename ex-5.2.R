rm(list=ls())
library(geoR)
data(SIC)
 points(sic.100, borders = sic.borders, pch=20,cex.max=3)


max.dist<-220
sic.bin<-variog(sic.100,option="bin",estimator.type="classical",max.dist=max.dist)
sic.rob<-variog(sic.100,option="bin",estimator.type="modulus",max.dist=max.dist)
plot(sic.bin$u,sic.bin$v,pch=20,col=1,ylab="semi-variogramme",
      xlab='h',ylim=c(0,max(c(sic.bin$v,sic.rob$v))))
points(sic.rob$u,sic.rob$v,pch=21,col=1)
legend(150, 5000, legend=c("classique","robuste"), col = c(1,1), pch=c(20,21))


plot(sic.bin$u,sic.bin$v,pch=20,ylab="semi-variogramme",
      xlab='h',ylim=c(0,max(c(sic.bin$v,sic.rob$v))))
ini<-c(15000,50) 
kappa<-0.5
cov.model<-"matern"
 wls <- variofit(sic.bin, ini = ini, cov.model=cov.model,weights="cressie",fix.nugget=TRUE,fix.kappa=FALSE,kappa=0.5)
 ols <- variofit(sic.bin, ini = ini, cov.model=cov.model,weights="equal",fix.nugget=TRUE,fix.kappa=FALSE,kappa=0.5)
 ml <- likfit(sic.100, ini = ini, cov.model=cov.model,fix.nugget=TRUE,fix.kappa=FALSE,kappa=0.5)
 lines(ols)
 lines(wls,lty=2)
  lines(ml,lty=3)
 legend(150, 5000, legend=c("MCO","MCP","MV"), lty=c(1,2,3))


set.seed(19)
wls.env <- variog.model.env(sic.100, obj.v = sic.bin,
                                  model.pars = wls,nsim=40)
     plot(sic.bin, env = wls.env,pch=20)


sic.xv.ols <- xvalid(sic.100, model = ols)
sic.xv.wls <- xvalid(sic.100, model = wls)
sic.xv.ml <- xvalid(sic.100, model = ml)
eqnm.ols<-mean(sic.xv.ols$std.error^2)
eqnm.wls<-mean(sic.xv.wls$std.error^2)
eqnm.ml<-mean(sic.xv.ml$std.error^2)


plot(sic.xv.wls, ask = FALSE,error = FALSE, std.error = TRUE, data.predicted = FALSE,pp = FALSE, map = FALSE, histogram = FALSE,error.predicted = TRUE,error.data = FALSE,pch=20)


ngridx<-100
ngridy<-100
xgrid<-seq(min(sic.borders[,1]),max(sic.borders[,1]),l=ngridx)
ygrid<-seq(min(sic.borders[,2]),max(sic.borders[,2]),l=ngridy)
pred.grid <- expand.grid(xgrid,ygrid)
krige.par<-krige.control(type.krige='ok',cov.pars=wls$cov.pars,cov.model=wls$cov.model)
kaq<-krige.conv(sic.100,locations=pred.grid,krige=krige.par)  
cols<- gray(6:16/16)
library(fields)
image.plot(xgrid,ygrid,matrix(kaq$predict,ngridx,ngridy),xlab='x',ylab='y',col =cols)
points(sic.100, pch=20,cex.max=1,pt.divide='equal',add.to.plot=TRUE, lty=2)
contour(kaq,add=T,nlevels=10)


image.plot(xgrid,ygrid,matrix(kaq$predict,ngridx,ngridy),xlab='x',ylab='y',col =cols)
points(sic.100, pch=20,cex.max=1,pt.divide='equal',add.to.plot=TRUE, lty=2)
contour(kaq,values=sqrt(kaq$krige.var),add=T,nlevels=10)


max.dist<-220
sic.cloud<-variog(sic.100,option="cloud",max.dist=max.dist)
sic.bin<-variog(sic.100,option="bin",max.dist=max.dist)
plot(sic.cloud,pch='.')
abline(v=sic.bin$bins.lim)
points(sic.bin$u,sic.bin$v,pch=20,cex=1.2)


