library(spatstat)
"newton.raphson"<-
function(xobs, xsim, psi, start = psi, maxiter = 30, eps1 = 1e-15, eps2 = 1e-08)
{

	iter <- 0
	nxt <- start - psi

	av <- apply(xsim, 2, mean)
	xobs <- xobs - av
	xsim <- sweep(xsim, 2, av)
	ll <- 0
	repeat {
		iter <- iter + 1
		cur <- nxt

		prob <- exp(xsim %*% cur)
		prob <- prob/sum(prob)

		E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)

		vtmp <- sweep(sweep(xsim, 2, E, "-"), 1, prob^0.5, "*")
		V <- t(vtmp) %*% vtmp

		nxt <- cur + (delt <- solve(V, xobs - E))
		ll.old <- ll
		repeat {
			ll <- sum(xobs * nxt) - log(mean(exp(xsim %*% nxt)))
					
			if(ll > ll.old - eps1)
				break
			else nxt <- cur + (delt <- delt/2)
		}
		if((abs(ll - ll.old) < eps1) || max(abs(delt/cur)) < eps2 || (
			iter >= maxiter))
			break
	}
	loglik <- ll.old

	cur <- nxt
        list(theta = cur + psi, se = diag(solve(V))^0.5, psi = psi, iter
			 = iter, loglik = loglik, E = E + av, V = V)
}

mcmcmle<-function(Q, interaction=NULL, covariates=NULL, correction="border", rbord=0, nsim=100, nrmh=1e5,
        start=list(n.start=X$n),control=list(nrep=nrmh),verb=TRUE,maxiter=30) {

#Q 	A data point pattern (of class "ppp") to which the model will be fitted, or a quadrature scheme (of class "quad") containing this pattern.
#interaction 	An object of class "interact" describing the point process interaction structure, or NULL indicating that a Poisson process (stationary or nonstationary) should be fitted.
#correction 	The name of the edge correction to be used. The default is "border" indicating the border correction. Other possibilities may include "Ripley", "isotropic", "translate" and "none", depending on the interaction.
#rbord 	If correction = "border" this argument specifies the distance by which the window should be eroded for the border correction.
#nsim 	Number of simulated realisations to generate 
#nrmh 	Number of Metropolis-Hastings iterations for each simulated realisation (for method="ho")
#start,control 	Arguments passed to rmh controlling the behaviour of the Metropolis-Hastings algorithm (for method="ho")
#verb 	Logical flag indicating whether to print progress reports (for method="ho")

fit.mpl<- ppm(Q,interaction=interaction, correction=correction, rbord=rbord, method="mpl")	
psi<- coef(fit.mpl)
k<- length(psi)
suff<-suffstat(fit.mpl, Q)
suffsim<-matrix(0,nsim,k)
    if (verb) 
        cat("Simulating... ")

for (i in 1:nsim) {
    if (verb) 
            cat(paste(i, " ", if (i%%10 == 0) 
                "\n", sep = ""))
    Xsim <- rmh(fit.mpl, control=list(nrep=nrmh),verbose=FALSE)
    suffsim[i,] <- suffstat(fit.mpl, Xsim)
    }    
    newton.raphson(suff,suffsim, psi, start = psi, maxiter = maxiter, eps1 = 1e-15, eps2 = 1e-08)
}
