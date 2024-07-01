#
# data available at http://www.image.ucar.edu/GSP/Data/O3.shtml
# Actual ozone measurements (max 8-hr daily average). Units are concentrations in parts per billion (PPB)
# The ozone "season" runs from 1st April through 1st October ( 184 days) over five years (1995 - 1999).

#
rm(list=ls())

source('ex-5.4-suppl.R')

ozmax8 <- matrix( scan("ozmax8.dat"), 920,513) 
#
# The file is sorted column-wise (first column is for station 1, etc). 
#
source("ozmax8.info.q") 


library(sp)
library(maps)

# we select a state

state<-'ohio'

a<- map('state', state,fill = TRUE,col=0)
title(paste(state,"state",sep=" "))
sel<-point.in.polygon(ozmax8.info$loc[,1], ozmax8.info$loc[,2], a$x, a$y)
sel <- (sel == 1)
state1<-'pennsylvania'

a<- map('state', state1,fill = TRUE,col=0)
title(paste(state1,"state1",sep=" "))
sel1<-point.in.polygon(ozmax8.info$loc[,1], ozmax8.info$loc[,2], a$x, a$y)
sel1 <- (sel1 == 1)
state<-c(state,state1)
sel <- (sel | sel1)

state2<-"west virginia"

a<- map('state', state2,fill = TRUE,col=0)
sel2<-point.in.polygon(ozmax8.info$loc[,1], ozmax8.info$loc[,2], a$x, a$y)
sel2 <- (sel2 == 1)
state<-c(state,state2)
sel <- (sel | sel2)


library(geoR)
a<-map('state', state,plot=FALSE)
coords<-ozmax8.info$loc[sel,]
obs<-ozmax8[,sel]
obs<-t(obs)
zs<-apply(obs,1,mean,na.rm=TRUE)
zt<-apply(obs,2,mean,na.rm=TRUE)
oz.geo<-as.geodata(cbind(coords,zs))
 points(oz.geo,xlab='longitude',ylab="latitude",cex.max=2)
lines(a$x,a$y)

# conversion

coords <- lonlat.to.planar(coords)



nsites<-dim(obs)[1]
nt<-dim(obs)[2]

z<-as.vector(obs)

tcoords<-(1:nt)%x%rep(1,dim(coords)[1])
xycoords<-(rep(1,nt)%x%as.matrix(coords))


# number of harmonics

ns<-6
h<-matrix(0,length(tcoords),2*ns)
tscale<-92
stations<-as.factor(rep(1,nt)%x%as.matrix(1:nsites))
for ( i in 1: ns) {
  h[,(2*(i-1)+1)]<-sin(pi*(i)*tcoords/tscale)
  h[,(2*i)]<-cos(pi*(i)*tcoords/tscale)
}
data.tmp<-data.frame(z=z,h=h,stations=stations)
fit.lm<-lm(z~h)
fitted<-predict(fit.lm,newdata=data.tmp)
znew<-z-fitted


site<-1
plot.ts(zt,ylab='PPB',xlab='Temps',lty=3)
lines(fitted[seq(site,length(z),nsites)])


site<-round(nsites/2,0)
plot.ts(z[seq(site,length(z),nsites)],ylab='PPB',xlab='Temps',lty=3)
title(main=paste("Station ", site, sep = ""))
lines(fitted[seq(site,length(z),nsites)])


data <-cbind(z,znew,xycoords,tcoords)
ntotal<-max(data[,5])
nsites<-dim(data)[1]/ntotal
coords<-data[1:nsites,3:4]




#
#
# we select the last years omitting October
#
sel<-(data[,5] > 736) & (data[,5] < 890)
z<-data[sel,2]

u <- as.vector(dist(as.matrix(coords)))
umax <- max(u)
tcoords<-data[sel,5]
xycoords<-data[sel,3:4]

lims<-define.bins(max.dist = umax,uvec=13)$bins.lim

tumax<-10+0.5
tlims<-c(0,seq(1.5,tumax,by=1))
tlims<-c(0,tlims)
sel<-!is.na(z)        
xycoords<-xycoords[sel,]
tcoords<-tcoords[sel]
znew<-z[sel]
binned<-empvar(znew,xycoords,tcoords,lims,tlims,maxdist=umax,tmaxdist=tumax)
variance <- var(znew)
scale <- c(umax,10)



cov.model <- "nsst"
q.nsst<-c(a=1.7,b=2,c=0.5,d=0.5,e=1,f=2)

lower <-c(0.00000001,0.000000001,0.0001,0.0001)
upper <-c(2*variance,2,2,0.999999)

theta<-c(variance,q.nsst[c(1,3,4)])
fit.scale <- TRUE
if (fit.scale) 
	{
	theta <- c(theta,scale)
	lower <-c(lower,0.000000001,0.000000001)
	upper <-c(upper,Inf,Inf)

	}
#
fit.nsst<-optim(theta,wls,lower = lower, upper = upper, binned=binned,cov.model= cov.model, param =q.nsst, scale = scale, fit.scale = fit.scale, weights='cressie', method= "L-BFGS-B",control=list(maxit=1000,parscale=theta))

# iacocesare 

cov.model <- "iacocesare"

q.iacocesare<-c(a=1.5,b=1.5,c=1.5,d=3)
theta<-c(variance,q.iacocesare[1:2])
lower <-c(0.000001,0,0)
upper <-c(2*variance,2,2)


fit.scale <- TRUE
if (fit.scale) 
	{
	theta <- c(theta,scale)
	lower <-c(lower,0.000000001,0.000000001)
	upper <-c(upper,Inf,Inf)

	}
#
fit.iacocesare<-optim(theta,wls,lower = lower, upper = upper, binned=binned,cov.model= cov.model, param =q.iacocesare, scale = scale, fit.scale = fit.scale, weights='cressie', method= "L-BFGS-B",control=list(maxit=1000,parscale=theta))
# gneiting14

cov.model <- "gneiting14"

q.gneiting14<-c(a=0.5,b=0.7,c=0.5,f=1.5,d=2)
theta<-c(variance,q.gneiting14[1:2])
lower <-c(0.000001,0.00000001,0)
upper <-c(4*variance,1,1)


fit.scale <- TRUE
if (fit.scale) 
	{
	theta <- c(theta,scale)
	lower <-c(lower,0.000001,0.0000001)
	upper <-c(upper,Inf,Inf)

	}
#

fit.gneiting14<-optim(theta,wls,lower = lower, upper = upper, binned=binned,cov.model= cov.model, param =q.gneiting14, scale = scale, fit.scale = fit.scale, weights='cressie', method= "L-BFGS-B",control=list(maxit=1000,parscale=theta))
# gneiting14sep

cov.model <- "gneiting14sep"
q.gneiting14sep<-c(a=0.5,b=0,c=0.5,f=1.5,d=2)
theta<-c(variance,q.gneiting14sep[1])
lower <-c(0.000001,0.00000001)
upper <-c(4*variance,1)

fit.scale <- TRUE
if (fit.scale) 
	{
	theta <- c(theta,scale)
	lower <-c(lower,0.000001,0.0000001)
	upper <-c(upper,Inf,Inf)

	}
#
fit.gneiting14sep<-optim(theta,wls,lower = lower, upper = upper, binned=binned,cov.model= cov.model, param =q.gneiting14sep, scale = scale, fit.scale = fit.scale, weights='cressie', method= "L-BFGS-B",control=list(maxit=1000,parscale=theta))

print(fit.nsst)
print(fit.iacocesare)
print(fit.gneiting14)
print(fit.gneiting14sep)




lag <- 0
plot(binned$centers,binned$emp.vario[,lag+1],xlab="distance", ylab="Semi-variogram",pch=20,ylim=c(0,1.3*variance))
h <- binned$centers
u<-binned$tcenters[lag+1]
title(main=paste("u = ",u,sep=""))

cov.model <- "nsst"
q<-c(a=NA,b=1,c=NA,d=NA,e=1,f=2)

param.nsst<-assign.param(fit.nsst$par, cov.model , param=q.nsst, fit.scale=TRUE)
semivar<- variog.theo(h,u,param.nsst)
lines(h,semivar,lty=4)


cov.model <- "iacocesare"
q<-c(a=NA,b=NA,c=1.5,d=3)
param.iacocesare<-assign.param(fit.iacocesare$par, cov.model , param=q, fit.scale=TRUE)
semivar<- variog.theo(h,u,param.iacocesare)
lines(h,semivar,col=1,lty=1)

cov.model <- "gneiting14"

q<-c(a=NA,b=NA,c=0.5,f=1.5,d=2)
param.gneiting14<-assign.param(fit.gneiting14$par, cov.model , param=q, fit.scale=TRUE)
semivar<- variog.theo(h,u,param.gneiting14)
lines(h,semivar,lty=2)

cov.model <- "gneiting14sep"
q<-c(a=NA,b=0,c=0.5,f=1.5,d=2)
param.gneiting14sep<-assign.param(fit.gneiting14sep$par, cov.model , param=q, fit.scale=TRUE)
semivar<- variog.theo(h,u,param.gneiting14sep)
lines(h,semivar,lty=3)
legend(600,100,legend=c("A","B","C","D"),lty=c(1,2,3,4))
print(param.nsst)
print(param.iacocesare)
print(param.gneiting14)
print(param.gneiting14sep)



plot(binned$tcenters,binned$emp.vario[1,],xlab="lag", ylab="Semi-variogram",pch=20,ylim=c(0,1.3*variance))
u <- binned$tcenters
h<-binned$centers[1]
title(main=paste("h = ",h,sep=""))

semivar<- variog.theo(h,u,param.nsst)
lines(u,semivar,lty=4)



semivar<- variog.theo(h,u,param.iacocesare)
lines(u,semivar,lty=1)


semivar<- variog.theo(h,u,param.gneiting14)
lines(u,semivar,lty=2)

semivar<- variog.theo(h,u,param.gneiting14sep)
lines(u,semivar,lty=3)
legend(7.5,100,legend=c("A","B","C","D"),lty=c(1,2,3,4))




lag <- 1
plot(binned$centers,binned$emp.vario[,lag+1],xlab="distance", ylab="Semi-variogram",pch=20,ylim=c(0,1.3*variance))
h <- binned$centers
u<-binned$tcenters[lag+1]
title(main=paste("u = ",u,sep=""))

semivar<- variog.theo(h,u,param.nsst)
lines(h,semivar,lty=4)


semivar<- variog.theo(h,u,param.iacocesare)
lines(h,semivar,col=1,lty=1)

semivar<- variog.theo(h,u,param.gneiting14)
lines(h,semivar,lty=2)

semivar<- variog.theo(h,u,param.gneiting14sep)
lines(h,semivar,lty=3)
legend(600,100,legend=c("A","B","C","D"),lty=c(1,2,3,4))




lag <- 2
plot(binned$centers,binned$emp.vario[,lag+1],xlab="distance", ylab="Semi-variogram",pch=20,ylim=c(0,1.3*variance))
h <- binned$centers
u<-binned$tcenters[lag+1]
title(main=paste("u = ",u,sep=""))

semivar<- variog.theo(h,u,param.nsst)
lines(h,semivar,lty=4)


semivar<- variog.theo(h,u,param.iacocesare)
lines(h,semivar,col=1,lty=1)

semivar<- variog.theo(h,u,param.gneiting14)
lines(h,semivar,lty=2)

semivar<- variog.theo(h,u,param.gneiting14sep)
lines(h,semivar,lty=3)

legend(600,100,legend=c("A","B","C","D"),lty=c(1,2,3,4))


# we use the last three days
sel<-(data[,5] > 0) & (data[,5] < 4)
sel.1<-(data[,5] == 4)



mse<-NULL
param<-param.nsst
mat<-stkrig.weights(data[sel.1,3:4],data[sel.1,5],data[sel,3:4],data[sel,5], param)
#
#
# we select the last three days of July
#
sse<-0
n<-0
for ( i in 1:31) {

sel<-(data[,5] > (885+i)) & (data[,5] < (889+i))
sel.1<-data[,5] == (889+i)
z.tmp<-data[sel,2]
missing<-is.na(z.tmp)
inv.cov<-solve(mat$covz[!missing,!missing])
lambda<-mat$c0[,!missing]%*%inv.cov
pred<-lambda%*%z.tmp[!missing]


error2<-(data[sel.1,2]-pred)^2
n<-n+sum(!is.na(error2))
sse<-sse+sum(error2,na.rm=TRUE)
}
mse<-c(mse,sse/n)

param<-param.iacocesare
mat<-stkrig.weights(data[sel.1,3:4],data[sel.1,5],data[sel,3:4],data[sel,5], param)
#
#
# we select the last three days on July
#
sse<-0
n<-0
for ( i in 1:31) {

sel<-(data[,5] > (885+i)) & (data[,5] < (889+i))
sel.1<-data[,5] == (889+i)
z.tmp<-data[sel,2]
missing<-is.na(z.tmp)
inv.cov<-solve(mat$covz[!missing,!missing])
lambda<-mat$c0[,!missing]%*%inv.cov
pred<-lambda%*%z.tmp[!missing]


error2<-(data[sel.1,2]-pred)^2
n<-n+sum(!is.na(error2))
sse<-sse+sum(error2,na.rm=TRUE)
}
mse<-c(mse,sse/n)

param<-param.gneiting14
mat<-stkrig.weights(data[sel.1,3:4],data[sel.1,5],data[sel,3:4],data[sel,5], param)
#
#
# we select the last three days on July
#
sse<-0
n<-0
for ( i in 1:31) {

sel<-(data[,5] > (885+i)) & (data[,5] < (889+i))
sel.1<-data[,5] == (889+i)
z.tmp<-data[sel,2]
missing<-is.na(z.tmp)
inv.cov<-solve(mat$covz[!missing,!missing])
lambda<-mat$c0[,!missing]%*%inv.cov
pred<-lambda%*%z.tmp[!missing]


error2<-(data[sel.1,2]-pred)^2
n<-n+sum(!is.na(error2))
sse<-sse+sum(error2,na.rm=TRUE)
}
mse<-c(mse,sse/n)

param<-param.gneiting14sep
mat<-stkrig.weights(data[sel.1,3:4],data[sel.1,5],data[sel,3:4],data[sel,5], param)
#
#
# we select the last three days on July
#
sse<-0
n<-0
for ( i in 1:31) {

sel<-(data[,5] > (885+i)) & (data[,5] < (889+i))
sel.1<-data[,5] == (889+i)
z.tmp<-data[sel,2]
missing<-is.na(z.tmp)
inv.cov<-solve(mat$covz[!missing,!missing])
lambda<-mat$c0[,!missing]%*%inv.cov
pred<-lambda%*%z.tmp[!missing]


error2<-(data[sel.1,2]-pred)^2
n<-n+sum(!is.na(error2))
sse<-sse+sum(error2,na.rm=TRUE)
}
mse<-c(mse,sse/n)


print(mse)


