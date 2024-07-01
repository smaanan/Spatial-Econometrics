rm(list=ls())
library(spdep)
data(eire)
lA <- lag.listw(nb2listw(eire.nb), eire.df$A)
set.seed(19)
a<-moran.test(spNamedVec("A", eire.df), nb2listw(eire.nb,style="W"),randomisation=FALSE,alternative="greater")


 colw <- nb2listw(eire.nb, ,style="W")
 nsim <- 999
 set.seed(1234)
 sim.M <- moran.mc(spNamedVec("A", eire.df), listw=colw, nsim=nsim,alternative="greater")

 b<-geary.test(spNamedVec("A", eire.df), nb2listw(eire.nb,style="W"),randomisation=FALSE,alternative="less")
 sim.G <- geary.mc(spNamedVec("A", eire.df), nb2listw(eire.nb, style="W"),
      nsim=nsim, alternative="less")
      


ra<- 27-rank(eire.df$towns)
A.lm <- lm(A ~ towns + pale, data=eire.df)
summary(A.lm)


res <- residuals(A.lm)
brks <- c(min(res),-2,-1,0,1,2,max(res))
         cols<- gray(0:5 / 5)
	 par(mfrow=c(1,1),pty='s')
moran.res<-  lm.morantest(A.lm, nb2listw(eire.nb,style="W"),alternative="greater")

A.SAR <- errorsarlm(A ~  pale, data=eire.df,nb2listw(eire.nb,style="W"), method="eigen", quiet=TRUE)
summary(A.SAR)
res <- residuals(A.SAR)
