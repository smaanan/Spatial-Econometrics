rm(list=ls())
y<-scan("coords-indigo.txt")
n<-649;x<-matrix(y[1:(3*n)],n,3,byrow=TRUE)

indigo<-read.csv('indigo.txt')
# we remove two outliers
indigo<-indigo[-c(380,381),]
indigo$y<-indigo$obs*2

ind<- sort(x[,1],index.return=TRUE)$ix
x<-x[ind,]
n<-dim(indigo)[1]
sel<- match(indigo$pt,x[,1])
indigo$xcoord<-x[sel,3]
indigo$ycoord<-x[sel,2]

sel<-c(56,300,451)

plot(indigo$xcoord[-sel],indigo$ycoord[-sel],pch='*',cex=0.5,xlab='longitude',ylab='latitude')
points(indigo$xcoord[sel],indigo$ycoord[sel],pch=20,cex=1.0)
text(indigo$xcoord[sel],indigo$ycoord[sel],pos=1)




test.indigo<-indigo[sel,]


library(mgcv)
# best model
indigo.gam<-gam(y~s(xcoord,ycoord)+GK+MK+TK+UK,data=indigo,family="poisson",subset=(-sel))
indigo.glm<-glm(y~GK+MK+TK+UK,data=indigo,family="poisson",subset=(-sel))

print(summary(indigo.gam))
print(summary(indigo.glm))
library(geoR)
library(geoRglm)
library(coda)
indigo.geo.res<-as.geodata(cbind(indigo$xcoord[-sel],indigo$ycoord[-sel],resid(indigo.glm)))
plot(variog(indigo.geo.res))



indigo.geo<-as.geodata(indigo[-sel,],coords.col=c(24,25),data.col=23,covar.col=c(8:11,20))
indigo.geo.test<-as.geodata(indigo[sel,],coords.col=c(24,25),data.col=23,covar.col=c(8:11,20))
covario <- covariog(indigo.geo,max.dist=0.4)  # sample covariogram
plot(covario)                      
parmval <- list(cov.model = "exponential", cov.pars = c(1, 0.05),max.dist =0.8)
class(parmval) <- "covariomodel"
lines(parmval, lty = 2)



trend.d<-trend.spatial("cte", indigo.geo, add=~GK+MK+TK+UK)
trend.l<-trend.spatial("cte", indigo.geo.test, add=~GK+MK+TK+UK)

locations <-  indigo[sel,24:25]
model<-model.glm.control(trend.d = trend.d, trend.l = trend.l, cov.model = "exponential")
set.seed(910)
output<-output.glm.control( sim.predict=TRUE, keep.mcmc.sim=TRUE)
# some tuning is necessary
a<-0.03
b<-0.21

phi.discrete <- seq(a, b, l =40)
prior<- prior.glm.control(phi.discrete = phi.discrete)
phi.scale<-phi.discrete[2]-phi.discrete[1]

# optimal scale proposal

S.scale <- 2*0.03
phi.scale  <-0.10*var(phi.discrete)
print(c(S.scale,phi.scale))
thin<-50
n.sample <-1000
 burn.in <- 1000
mcmc.input <- mcmc.control(S.scale = S.scale, burn.in = burn.in, n.iter = n.sample*thin, thin = thin, phi.scale = phi.scale, phi.start=mean(phi.discrete))

fit.bayes <- pois.krige.bayes(indigo.geo, locations = locations, prior = prior, mcmc.input = mcmc.input, model = model,, output =output)
# diagnostics
library(coda)
indigo.coda<-create.mcmc.coda(fit.bayes, mcmc.input = mcmc.input)

para<-fit.bayes$posterior$beta$sample[1,]
d<-density(para)
h<-hist(para,freq=FALSE,plot=FALSE)
ylim<-c(0,max(max(h$den),max(d$y)))
hist(para,freq=FALSE,main=' ',xlab = 'cte',ylim=ylim)
lines(d)



para<-fit.bayes$posterior$beta$sample[2,]
d<-density(para)
h<-hist(para,freq=FALSE,plot=FALSE)
ylim<-c(0,max(max(h$den),max(d$y)))
hist(para,freq=FALSE,main=' ',xlab = 'GK',ylim=ylim)
lines(d)



para<-fit.bayes$posterior$beta$sample[3,]
d<-density(para)
h<-hist(para,freq=FALSE,plot=FALSE)
ylim<-c(0,max(max(h$den),max(d$y)))
hist(para,freq=FALSE,main=' ',xlab = 'MK',ylim=ylim)
lines(d)



para<-fit.bayes$posterior$beta$sample[4,]
d<-density(para)
h<-hist(para,freq=FALSE,plot=FALSE)
ylim<-c(0,max(max(h$den),max(d$y)))
hist(para,freq=FALSE,main=' ',xlab = 'TK',ylim=ylim)
lines(d)



para<-fit.bayes$posterior$beta$sample[5,]
d<-density(para)
h<-hist(para,freq=FALSE,plot=FALSE)
ylim<-c(0,max(max(h$den),max(d$y)))
hist(para,freq=FALSE,main=' ',xlab = 'UK',ylim=ylim)
lines(d)



para<-fit.bayes$posterior$sigma$sample
d<-density(para)
h<-hist(para,freq=FALSE,plot=FALSE)
ylim<-c(0,max(max(h$den),max(d$y)))
hist(para,freq=FALSE,main=' ',xlab = expression(sigma^2),ylim=ylim)
lines(d)




para<-fit.bayes$posterior$phi$sample
d<-density(para)
h<-hist(para,freq=FALSE,plot=FALSE)
ylim<-c(0,max(max(h$den),max(d$y)))
hist(para,freq=FALSE,main=' ',xlab = expression(phi),ylim=ylim)
lines(d)




lambda<-fit.bayes$predictive$simulations[1,]
pred<-rpois(length(lambda),lambda)
plot(table(pred)/length(pred),xlab='site 1',ylab='prob')
points(indigo$y[sel][1],0,pch=20,cex=4)



lambda<-fit.bayes$predictive$simulations[2,]
pred<-rpois(length(lambda),lambda)
plot(table(pred)/length(pred),xlab='site 2',ylab='prob')
points(indigo$y[sel][2],0,pch=20,cex=4)



lambda<-fit.bayes$predictive$simulations[3,]
pred<-rpois(length(lambda),lambda)
plot(table(pred)/length(pred),xlab='site 3',ylab='prob')
points(indigo$y[sel][3],0,pch=20,cex=4)


