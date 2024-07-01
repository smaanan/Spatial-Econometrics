
pcp.sim.new<-function (rho, m, s2, region.poly, larger.region = NULL, vectorise.loop = TRUE) 
{
    if (is.null(larger.region)) 
        larger.region <- as.vector(apply(sbox(region.poly), 2, 
            range))
    sim.events <- c(0, 0)
    n <- rpois(1, lambda = rho * (larger.region[2] - larger.region[1]) * 
        (larger.region[4] - larger.region[3]))
    parents <- cbind(runif(n, larger.region[1], larger.region[2]), 
        runif(n, larger.region[3], larger.region[4]))
    sd <- sqrt(s2)
    if (!vectorise.loop) {
        for (j in 1:n) {
            num.child <- rpois(1, lambda = m)
            for (k in 1:num.child) {
                new.child <- parents[j, ] + rnorm(2, 0, sd = sd)
                sim.events <- rbind(sim.events, new.child)
            }
        }
        sim.events <- sim.events[-1, ]
    }
    else {
        num.child <- rpois(n, lambda = m)
        num.tot <- sum(num.child)
        sim.events <- matrix(c(rnorm(num.tot, 0, sd) + rep(parents[, 
            1], num.child), rnorm(num.tot, 0, sd) + rep(parents[, 
            2], num.child)), num.tot)
    }
    child<-sim.events
    return(list(child=child,parents=parents))
}


# plot the data

library(spatstat)
library(splancs)

data(finpines)

X<-finpines
x<-X$x
y<-X$y
plot(x,y,xlab="x",ylab="y",pch=20)


# CSR envelope 

poly<-rbind(c(X$w$xrange[1],X$w$yrange[1]),c(X$w$xrange[2],X$w$yrange[1]),c(X$w$xrange[2],X$w$yrange[2]),c(X$w$xrange[1],X$w$yrange[2]))
nsim<-40
a<-Kest(X)
r.range<-attributes(a)$alim
r <- seq(r.range[1],r.range[2],length=180)
K.env <- Kenv.csr(length(x), poly, nsim=nsim, r)
L.env <- lapply(K.env, FUN=function(x) sqrt(x/pi)-r)
K.est<-khat(as.points(X), poly, r)
L.est<- sqrt(K.est/pi)-r
limits <- range(c(unlist(L.env),L.est))

plot(r, L.est, type="l", xlab="h",ylim=limits,ylab=expression(hat(L)(h)-h))
lines(r, L.env$upper, lty=2)
lines(r, L.env$lower, lty=2)
abline(h=0)
# fit a model
pcp.fit <- pcp(cbind(x,y), poly, h0=r.range[2], n.int=100)
m <- npts(cbind(x,y))/(areapl(poly)*pcp.fit$par[2])
# simulate a realization form the fitted model
sims <- pcp.sim.new(pcp.fit$par[2], m, pcp.fit$par[1], poly)


plot(sims$child[,1],sims$child[,2],xlab="x",ylab="y",pch=20,cex=0.6)
points(sims$parents[,1],sims$parents[,2],pch=2,cex=1.0)


K.env <- Kenv.pcp(pcp.fit$par[2], m, pcp.fit$par[1], poly, nsim=nsim, r=r)
L.env <- lapply(K.env, FUN=function(x) sqrt(x/pi)-r)
limits <- range(c(unlist(L.env),L.est))
plot(r, L.est, ylim=limits,type="l",
     xlab="h",ylab=expression(hat(L)(h)-h))
#define the theoretical K function for a Neyman-Scott process
  K.NS<-function(h,s2,ro)
    {
   return(pi*h*h+(1-exp(-h*h*0.25/s2))/ro)
    }
     
K.theo<-K.NS(r,pcp.fit$par[1],pcp.fit$par[2])
L.theo<-sqrt(K.theo/pi)-r     
lines(r, L.env$lower, lty=5)
lines(r, L.env$upper, lty=5)
lines(r, L.theo, lty=3,lwd=3)
abline(h=0)


