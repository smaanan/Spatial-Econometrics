library(geoR) #loads geoR into your R session




data(parana) #loads the data parana into your R session.
par(mar = c(3,3,1,1))
plot(variog(parana, option='cloud'),pch='.',cex=0.5)
a<-variog(parana, option='smooth', band=100)
lines(a$u,a$v)




par(mar = c(3,3,1,1))
a1<-variog(parana, option='smooth', direction=0,band=100)
a2<-variog(parana, option='smooth', direction=pi/4,band=100)
a3<-variog(parana, option='smooth', direction=pi/2,band=100)
a4<-variog(parana, option='smooth', direction=3*pi/4,band=100)
n.o<-c("0","45","90","135")
ymax <- max(c(a1$v, a2$v, a3$v, a4$v), na.rm = TRUE)

plot(a1$u,a1$v,type='l',xlab='distance',ylab='semivariance',ylim = c(0, ymax))

lines(a2$u,a2$v,lty=2)
lines(a3$u,a3$v,lty=3,lwd=2)

lines(a4$u,a4$v,lty=4)
legend(0, ymax, legend = c(substitute(a * degree, 
                  list(a = n.o[1])), substitute(a * degree, list(a = n.o[2])), 
                  substitute(a * degree, list(a = n.o[3])), substitute(a * 
                    degree, list(a = n.o[4])), expression()),lty = c(1:4),lwd=c(1,1,2,1))



par(mar = c(3,3,1,1),pch=20)
a1<-variog(parana, option='smooth', direction=0,band=100,trend='1st')
a2<-variog(parana, option='smooth', direction=pi/4,band=100,trend='1st')
a3<-variog(parana, option='smooth', direction=pi/2,band=100,trend='1st')
a4<-variog(parana, option='smooth', direction=3*pi/4,band=100,trend='1st')
n.o<-c("0","45","90","135")
ymax <- max(c(a1$v, a2$v, a3$v, a4$v), na.rm = TRUE)

plot(a1,type='l',xlab='distance',ylab='semivariance',ylim = c(0, ymax))
lines(a2$u,a2$v,lty=2)
lines(a3$u,a3$v,lty=3,lwd=2)

lines(a4$u,a4$v,lty=4)
legend(0, ymax, legend = c(substitute(a * degree, 
                  list(a = n.o[1])), substitute(a * degree, list(a = n.o[2])), 
                  substitute(a * degree, list(a = n.o[3])), substitute(a * 
                    degree, list(a = n.o[4])), expression()),lty = c(1:4),lwd=c(1,1,2,1))


par(mar = c(3,3,1,1),pch=20)
a1<-variog(parana, option='smooth', direction=0,band=100,trend='2nd')
a2<-variog(parana, option='smooth', direction=pi/4,band=100,trend='2nd')
a3<-variog(parana, option='smooth', direction=pi/2,band=100,trend='2nd')
a4<-variog(parana, option='smooth', direction=3*pi/4,band=100,trend='2nd')
n.o<-c("0","45","90","135")
ymax <- max(c(a1$v, a2$v, a3$v, a4$v), na.rm = TRUE)

plot(a1,type='l',xlab='distance',ylab='semivariance',ylim = c(0, ymax))
lines(a2$u,a2$v,lty=2)
lines(a3$u,a3$v,lty=3,lwd=2)

lines(a4$u,a4$v,lty=4)
legend(0, ymax, legend = c(substitute(a * degree, 
                  list(a = n.o[1])), substitute(a * degree, list(a = n.o[2])), 
                  substitute(a * degree, list(a = n.o[3])), substitute(a * 
                    degree, list(a = n.o[4])), expression()),lty = c(1:4),lwd=c(1,1,2,1))


par(mfrow=c(1,1),mar = c(3,3,1,1))
parana.cloud<-variog(parana, option='cloud')
parana.bin<-variog(parana,option="bin",estimator.type="classical",trend='cte')
plot(parana.cloud,ylab="semi-variogramme",xlab='h',pch='.',cex=0.3)
points(parana.bin$u,parana.bin$v,pch=20,col=1,cex=1.4)
abline(v=parana.bin$bins.lim)
text(parana.bin$uvec,parana.bin$v,label=parana.bin$n,pos=3,offset=1,cex=1.1)


par(mfrow=c(1,1),mar = c(3,3,1,1))
parana.bin<-variog(parana,option="bin",estimator.type="classical",trend='cte')
parana.rob<-variog(parana,option="bin",estimator.type="modulus",trend='cte')
plot(parana.bin$u,parana.bin$v,pch=20,col=1,ylab="semi-variogramme",
      xlab='h',ylim=c(0,max(c(parana.bin$v,parana.rob$v))))
points(parana.rob$u,parana.rob$v,pch=21,col=1)
legend(70, 6000, legend=c("classique","robuste"), col = c(1,1), pch=c(20,21))


par(mfrow=c(1,1),mar = c(3,3,1,1))
parana.bin<-variog(parana,option="bin",estimator.type="classical",trend='1st')
plot(parana.bin$u,parana.bin$v,pch=20,ylab="semi-variogramme",xlab='h',ylim=c(0,max(parana.bin$v)))
ini<-c(1000,70) 
cov.model<-"exp"
 wls <- variofit(parana.bin, ini = ini, cov.model=cov.model,weights="cressie")
 ols <- variofit(parana.bin, ini = ini, cov.model=cov.model,weights="equal")
 ml <- likfit(parana, method="ML",trend="1st", ini = ini, cov.model=cov.model)
 lines(ols)
 lines(wls,lty=2)
  lines(ml,lty=3)
 legend(400, 400, legend=c("MCO","MCP","MV"), lty=c(1,2,3))


cov.model<-"exp"
fix.kappa<-TRUE
kappa<-0.5
ml1 <- likfit(parana,cov.model=cov.model,nugget=500,ini=c(1000,70),method="ML",trend="1st",fix.kappa=fix.kappa,kappa=kappa)

ml1.aniso <- likfit(parana,cov.model=cov.model,nugget=500,ini=c(1000,70),method="ML",trend="1st", fix.psiA = FALSE, fix.psiR = FALSE,fix.kappa=fix.kappa,kappa=kappa)

ml2 <- likfit(parana,cov.model=cov.model,nugget=500,ini=c(100,20),method="ML",trend="2nd",fix.kappa=fix.kappa,kappa=kappa)

xv.ml1 <- xvalid(parana, model = ml1)
xv.ml1.aniso <- xvalid(parana, model = ml1.aniso)
xv.ml2 <- xvalid(parana, model = ml2)
eqnm.ml1<-mean(xv.ml1$std.error^2)
eqnm.ml1.aniso<-mean(xv.ml1.aniso$std.error^2)
eqnm.ml2<-mean(xv.ml2$std.error^2)


par(mar = c(3,3,0,1))
par(mgp = c(2,1,0))
plot(xv.ml2, ask = FALSE,error = FALSE, std.error = TRUE, data.predicted = FALSE,pp = FALSE, map = TRUE, histogram = FALSE,error.predicted = FALSE,error.data = FALSE,pch=20)


par(mar = c(3,3,0,1))
par(mgp = c(2,1,0))
plot(xv.ml2, ask = FALSE,error = FALSE, std.error = TRUE, data.predicted = FALSE,pp = FALSE, map = FALSE, histogram = FALSE,error.predicted = TRUE,error.data = FALSE,pch=20)



par(mar = c(3,3,1,1),pch=20)
nsim<-40
parana.vario <- variog(parana, trend='2nd')
parana.env <- variog.model.env(parana, obj.v = parana.vario,model.pars = ml2,nsim=nsim)
parana.vario <- variog(parana)
plot(parana.vario, env = parana.env)


par(mar = c(3,3,1,1),pch=20)
ngridx<-100
ngridy<-100
xgrid<-seq(min(parana$borders[,1]),max(parana$borders[,1]),l=ngridx)
ygrid<-seq(min(parana$borders[,2]),max(parana$borders[,2]),l=ngridy)
cols<-gray(seq(1,0.1,l=30))

#Universal Kriging
pred.grid <- expand.grid(x=xgrid,y=ygrid) # defines a grid, with (min,max,length of vector) of x coordinate, and (min, max, length of vector) of y coordinate.
krige.ml <- krige.conv(parana,loc=pred.grid,krige=krige.control(obj.m=ml2,type.krige='ok',trend.d="2nd",trend.l="2nd"),borders=NULL)
image(krige.ml,loc=pred.grid,col=gray(seq(1,0.1,l=30)), xlab="O-E", ylab="S-N",bor=parana$borders)
contour(krige.ml,loc=pred.grid,values=krige.ml$pred, xlab="O-E", ylab="S-N",bor=parana$borders,add=TRUE) 
points(parana$coords,pch=20,cex=0.7)





par(mar = c(3,3,1,1),pch=20)
image(krige.ml,loc=pred.grid,col=gray(seq(1,0.1,l=30)), xlab="O-E", ylab="S-N",bor=parana$borders) 
contour(krige.ml,loc=pred.grid,values=sqrt(krige.ml$krige.var), xlab="O-E", ylab="S-N",bor=parana$borders,add=TRUE,nlevels=15)
points(parana$coords,pch=20,cex=0.7)


     par(mar = c(3,3,1,1),pch=20)
library(fields)
ngridx<-100
ngridy<-100
xgrid<-seq(min(parana$borders[,1]),max(parana$borders[,1]),l=ngridx)
ygrid<-seq(min(parana$borders[,2]),max(parana$borders[,2]),l=ngridy)
cols<-gray(seq(1,0.1,l=30))
#Universal Kriging
pred.grid <- expand.grid(x=xgrid,y=ygrid) # defines a grid, with (min,max,length of vector) of x coordinate, and (min, max, length of vector) of y coordinate.

krige.ml <- krige.conv(parana,loc=pred.grid,krige=krige.control(obj.m=ml2,type.krige='ok',trend.d="2nd",trend.l="2nd"),borders=NULL)
image.plot(xgrid,ygrid,matrix(krige.ml$predict,ngridx,ngridy),xlab='x',ylab='y',col =cols)
points(parana, pch=20,cex.max=1,pt.divide='equal',add.to.plot=TRUE, lty=2)
contour(krige.ml,add=T,nlevels=10)


par(mar = c(3,3,1,1),pch=20)
image.plot(xgrid,ygrid,matrix(krige.ml$predict,ngridx,ngridy),xlab='x',ylab='y',col =cols)
points(parana, pch=20,cex.max=1,pt.divide='equal',add.to.plot=TRUE, lty=2)
contour(krige.ml,values=sqrt(krige.ml$krige.var),add=T,nlevels=10)


