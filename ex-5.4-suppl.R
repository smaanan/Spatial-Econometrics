dyn.load(paste("ex-5.4-suppl",.Platform$dynlib.ext,sep=""))
# Montserrat Fuentes' program  (http://www4.stat.ncsu.edu/~fuentes/)
# to compute geodetic distance

rdistearth<-
function(loc1, loc2, miles = FALSE )
{
	if (miles) 
            R <- 3963.34
        else R <- 6378.388	
        if(missing(loc2))
                loc2 <- loc1
        R <- 6371
        lat <- loc1[, 2]
        lon <- loc1[, 1]
        coslat1 <- cos((lat * pi)/180)
        sinlat1 <- sin((lat * pi)/180)
        coslon1 <- cos((lon * pi)/180)
        sinlon1 <- sin((lon * pi)/180)
        lat <- loc2[, 2]
        lon <- loc2[, 1]
        coslat2 <- cos((lat * pi)/180)
        sinlat2 <- sin((lat * pi)/180)
        coslon2 <- cos((lon * pi)/180)
        sinlon2 <- sin((lon * pi)/180)
        PP1 <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1)
        PP2 <- cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2)
        pp <- (PP1 %*% t(PP2))
        R * acos(ifelse(pp > 1, 1, pp))
}


# transform the coordinates
# Montserrat Fuentes' projection on the plane Centered around the center of gravity

lonlat.to.planar<-
function(lon.lat, miles =FALSE)
{
        x <- lon.lat[, 1]
        y <- lon.lat[, 2]
        mx <- mean(x)
        my <- mean(y)
        temp <- cbind(rep(mx, 2), range(y))
        sy <- rdistearth(temp)[2, 1]
        temp <- cbind(range(x), rep(my, 2))
        sx <- rdistearth(temp ,miles = miles)[2, 1]
        temp <- list(x = sx/(max(x) - min(x)), y = sy/(max(y) - min(y)))
        cbind((x - mx) * temp$x, (y - my) * temp$y)
}

define.bins<-function (max.dist, uvec = "default", breaks = "default", nugget.tolerance=0) 
{
    if (all(breaks == "default")) {
        if (all(uvec == "default")) 
            uvec <- 13
        if (mode(uvec) == "numeric") {
            if (length(uvec) == 1) {
                bins.lim <- seq(0, max.dist, l = uvec + 1)
                bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim > 
                  nugget.tolerance])
                uvec <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])
            }
            else {
                uvec <- c(0, uvec)
                nvec <- length(uvec)
                d <- 0.5 * diff(uvec[2:nvec])
                bins.lim <- c(0, (uvec[2:(nvec - 1)] + d), (d[nvec - 
                  2] + uvec[nvec]))
                bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim > 
                  nugget.tolerance])
            }
        }
        else stop("argument uvec can only take a numeric vector")
    }
    else {
        if (mode(breaks) != "numeric") 
            stop("argument breaks can only take a numeric vector")
        else bins.lim <- breaks
        bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim > 
            nugget.tolerance])
        uvec <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])
    }
    return(list(uvec = uvec, bins.lim = bins.lim))
}

# A large part of the following code has been adapted from the R package RandomFields
# author Martin Schlather http://www.stochastik.math.uni-goettingen.de/~schlather/
#
# Calculates the empirical (semi-)variogram of a spatio-temporal random field realisation
#
empvar<-function(z, coords, tcoords, lims, tlims, maxdist, tmaxdist,modulus =FALSE, sdcalc =TRUE) {
   if (lims[1] < 1e-16) 
      lims[1] <- -1
         if (tlims[1] < 1e-16) 
            tlims[1] <- -1
        n.data<-length(z)
        nbins<-length(lims)-1
        ntbins<-length(tlims)-1
	n<-(nbins)*(ntbins)
        vbin<-numeric(n)
        cbin<-numeric(n)
	sdbin<-numeric(n)
        
        a <- .C("empvar", n.data=as.integer(n.data), as.double(as.vector(coords[,1])), as.double(as.vector(coords[, 2])), 
                tc=as.double(as.vector(tcoords)), as.double(z), nbins = as.integer(nbins),ntbins = as.integer(ntbins), lims=as.double(as.vector(lims)), tlims=as.double(as.vector(tlims)), 
                as.integer(modulus), as.double(maxdist), as.double(tmaxdist), 
                cbin = as.integer(cbin), vbin = as.double(vbin), as.integer(sdcalc),sdbin = as.double(sdbin))

	 a$vbin[a$cbin == 0] <-NA
         a$vbin[1]<-0
         vv<-matrix(a$vbin,a$nbins,a$ntbins)
         cv<-matrix(a$cbin,a$nbins,a$ntbins)

# strange correction
# why ? I made this in order to get results consistent with Schlather.
	 cv[,1]<-2*cv[,1]


         sv<-matrix(a$sdbin,a$nbins,a$ntbins)
	bins<-lims[-1]

	tbins<-tlims[-1]
	tbins[1]<-0.5	
	centers<-bins-c(0,diff(bins)/2)
        tcenters<- tbins-c(0,diff(tbins)/2)
	tcenters[1]<-0
	return(list(emp.vario=vv,cv=cv,sv=sv,bins=bins, tbins=tbins,centers=centers,tcenters=tcenters))
        }

logit<-function(x) log(x/(1-x))
inv.logit<-function(x) exp(x)/(1+exp(x))


matern.fun<-function (u, phi, kappa) 
#
# geoR code
#
{
    if (is.vector(u)) 
        names(u) <- NULL
    if (is.matrix(u)) 
        dimnames(u) <- list(NULL, NULL)
    uphi <- u/phi
    uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/gamma(kappa)) * 
        (uphi^kappa) * besselK(x = uphi, nu = kappa)), 1)
    uphi[u > 600 * phi] <- 0
    return(uphi)
}

gneiting14.cor<-function(h,u,param) {
	
# eq. (14) in Gneting (2002)
# alpha = param[1]
# beta = param[2]
# gamma = param[3]
# tau = param[4]
# d = param[5] 
# The parameter param[5] is the genuine spatial dimension of the field
# check the parameters

# Changes in the future
#  2*param[1] -> param[1]
#  2*param[3] -> param[3]
#  2*param[3] -> param[3]
#  param[2] -> param[2]/2

	if ( (param[1] <=0) | (param[1] >1))
 		stop("invalid param[1] parameter in gneiting.cor function")
	if ( (param[2] <0) | (param[2] >1))
		stop("invalid param[2] parameter in gneiting.cor function")
	if ( (param[3] <=0) | (param[3] >1))
 		stop("invalid param[3] parameter in gneiting.cor function")
	if  (param[4] < param[5]/2)
 		stop("invalid param[4] or param[5] parameter in gneiting.cor function")
        b <-(1+u^(2*param[1]))

	cv<-exp(-(h^2/b^param[2])^param[3])/(b^param[4])

	return(cv)
}

iacocesare.cor<-function(h,u,param)
{
# check the parameter
	if ( (param[1] < 0) | (param[1] > 2))
 		stop("invalid param[1] parameter in iacocesare.cor function")
	if ( (param[2] < 0) | (param[2] > 2))
		stop("invalid param[2] parameter in iacocesare.cor function")
	if  (param[3] < param[4]/2)
 		stop("invalid param[3] or param[4] parameter in iacocesare.cor function")

	correlation <- (1+h^param[1]+u^param[2])^{-param[3]}
        return(correlation)

}


nsst.cor<-function(h,u,param)
# nsst (Non-Separable Space-Time model)

# C(x,t)= psi(t)^{-f} phi(x / psi(t))

#This model is used for space-time modelling where the spatial component is isotropic.
#phi is the stable model if b=1;
#phi is the whittlematern model if b=2;
#phi is the cauchy model if b=3;
#Here, a is the respective parameter for the model; the restrictions on a are described there.

#The function psi satisfies
#psi^2(t) = (t^c+1)^d if e=1
#psi^2(t) = (d^{-1} t^c+1)/(t^c+1) if e=2
#psi^2(t) = -log(t^c+1/d)/log d if e=3
#The parameter f must be greater than or equal to the genuine spatial dimension of the field. Furthermore, c in (0,2] and d in (0,1). The spatial dimension must be >=1. 

{

# Warning !!!!!
# The parameter f must be greater than or equal to the genuine spatial dimension of the field. Here we put this dimension
# equal to 2
#
	if ((param[6] < 2) )
		 		stop("invalid param[6] parameter in nsst.cor function")
        if (sum((param[5] == c(1,2,3))) <= 0)
					stop("invalid param[5] parameter in nsst.cor function")
        if (sum((param[2] == c(1,2,3))) <= 0)
					stop("invalid param[2] parameter in nsst.cor function")
	if ( (param[3] <= 0) | (param[3] > 2))
		stop("invalid param[1] parameter in nsst.cor function")
	if ( (param[4] <= 0) | (param[4] >= 1))
	 		stop("invalid param[4] parameter in nsst.cor function")
        
	if (param[5] == 1) {
		psi<-function(u,param)
		 { 
		 val <-(u^param[3]+1)^param[4]
                 return(sqrt(val))
		 }	
	}

       if (param[5] == 2) {
		psi<-function(u,param)
		 { 
		 val<-((u^param[3])/param[4]+1)/(u^param[3]+1)

                 return(sqrt(val))
		 }	
	}

       if (param[5] == 3) {
		psi<-function(u,param)
		 { 
		 val<--log(u^param[3]+1/param[4])/log(param[4])                      
                 return(sqrt(val))
		 }	
	}




       if (param[2] == 1) {
		if ( (param[1] <= 0) | (param[1] > 2))
	 		stop("invalid param[1] parameter in nsst.cor function")

		phi<-function(h,param)
		 { 
		 val <- exp(-h^param[1])
                 return(val)
		 }	
	}
       if (param[2] == 2) {
		if ( (param[1] <= 0))
	 		stop("invalid param[1] parameter in nsst.cor function")
		phi<-function(h,param)
		 { 

		 val <- matern.fun(h,1, param[1])
                 return(val)
		 }	
	}


       if (param[2] == 3) {
		if ( (param[1] <= 0))
	 		stop("invalid param[1] parameter in nsst.cor function")

		phi<-function(h,param)
		 { 
		  val <- (1+h^2)^(-param[1])
                 return(val)
		 }	
	}


	th <- h/psi(u,param)
        correlation <-  phi(th,param)/(psi(u,param)^param[6])
        return(correlation)
}





correlation.fun<-function(h,u,param,cov.model = NULL,scale = NULL)
{
       if (missing(scale)){ 
            scale <- c(1,1)
        }
        hs<-h/scale[1]
	us<-u/scale[2]

	cov.model <- match.arg(cov.model, choices = c("gneiting14","gneiting14sep","iacocesare", "nsst"))
	if (cov.model == "gneiting14sep") 
		
        	return(gneiting14.cor(hs,us,param))

	if (cov.model == "gneiting14") 
		
        	return(gneiting14.cor(hs,us,param))
        
	if (cov.model == "iacocesare") 

        	return(iacocesare.cor(hs,us,param))

	if (cov.model == "nsst") 

        	return(nsst.cor(hs,us,param))

}

cov.fun<-function(h,u,param,variance=1,cov.model = NULL,scale = NULL)
{
	val<-variance*correlation.fun(h,u,param,cov.model,scale)
	return(val)
}


variog.fun<-function(h,u,param,variance=1,cov.model = NULL,scale = NULL)
{
	val<-variance*(1-correlation.fun(h,u,param,cov.model,scale))
	return(val)
}

cov.theo<-function(h,u,a)
{
        variance <- a$variance
        param <- a$param
        scale <- a$scale 
	cov.model <-a$cov.model
	val <- cov.fun(h,u,param,variance,cov.model, scale)
	return(val)

}

variog.theo<-function(h,u,a)
{
        variance <- a$variance
        param <- a$param
        scale <- a$scale 
	cov.model <-a$cov.model
	val <- variog.fun(h,u,param,variance,cov.model, scale)
	return(val)

}

variog.fun.old<-function(h,u,param,variance=1,cov.model = NULL,scale = NULL)
{
        if (missing(scale)){ 
            scale <- c(1,1)
        }
        hs<-h/scale[1]
	us<-u/scale[2]

	cov.model <- match.arg(cov.model, choices = c("gneiting14","gneiting14sep","iacocesare", "nsst"))
	if (cov.model == "gneiting14sep") 
		
        	return(variance*(1-gneiting14.cor(hs,us,param)))

	if (cov.model == "gneiting14") 
		
        	return(variance*(1-gneiting14.cor(hs,us,param)))
        
	if (cov.model == "iacocesare") 

        	return(variance*(1-iacocesare.cor(hs,us,param)))

	if (cov.model == "nsst") 

        	return(variance*(1-nsst.cor(hs,us,param)))

}
assign.param <- function(theta, cov.model , param, fit.scale) {
  	if (cov.model == "gneiting14")
        	{
		variance<-theta[1]
		param[1:2]<-theta[2:3]
		if (fit.scale) 
	            scale<-theta[4:5]

		}
	if (cov.model == "gneiting14sep")
        	{
		variance<-theta[1]
		param[1]<-theta[2]
		if (fit.scale) 
	            scale<-theta[3:4]

		}

        if (cov.model == "iacocesare")
        {
		variance<-theta[1]
		param[1:2]<-theta[2:3]
		if (fit.scale) 
	            scale<-theta[4:5]

	}
       if (cov.model == "nsst")
        {
		variance<-theta[1]
		param[1]<-theta[2]
		param[3]<-theta[3]
		param[4]<-theta[4]
		if (fit.scale) 
	            scale<-theta[5:6]

	}
	
       return(list(variance=variance,param=param, scale=scale, cov.model =cov.model))	
}

wls<-function(theta, binned, cov.model = NULL, param = NULL, scale = NULL, fit.scale = FALSE, weights) {
        if (missing(weights)){ 
            weights <- "cressie"
        }
        else weights <- match.arg(weights, choices = c("npairs", 
        "equal", "cressie"))
        if (missing(scale)){ 
            scale <- c(1,1)
        }
        variog.param<-assign.param(theta, cov.model , param, fit.scale)
	loss<-0
	for ( i in 1:dim(binned$emp.vario)[2]) {
	sel<-(!is.na(binned$emp.vario[,i]))
        if ( i == 1 ) 
	    sel[1]<-FALSE
        gamma.hat<-as.numeric(binned$emp.vario[sel,i])
        u<-binned$tcenters[i]
        n<-binned$cv[sel,i]
        h<-as.numeric(binned$centers[sel])
	gamma.th <- variog.theo(h,u,variog.param)
	if (weights == "equal") 
             tloss <- sum((gamma.hat - gamma.th)^2)
        if (weights == "npairs") 
             tloss <- sum(n * (gamma.hat - gamma.th)^2)
        if (weights == "cressie") 
		tloss<-sum(n*(1-gamma.hat/gamma.th)^2)

        loss <- loss + tloss
	}

	return(loss)

}

stkrig.weights<-function(floc,ftimes,loc,times,  model)
#
# Computes the weights assign for each data point in simple  krigring
#
{
	variance <- model$variance
        param <- model$param
        scale <- model$scale 
	cov.model <-model$cov.model
	cov.param<-assign.param(theta, cov.model , param, fit.scale)	

        cov.param<-model 
	n<-dim(loc)[1]
	nn<-dim(floc)[1]	
        
	covz<-matrix(0,n,n)
	for ( i in 1: n) 
	     for ( j in i:n) {
	        	sdist<-sqrt((loc[i,1]-loc[j,1])^2+(loc[i,2]-loc[j,2])^2)           
			tdist<-abs(times[i]-times[j])
       covz[j,i]<-covz[i,j]<-cov.theo(sdist,tdist,cov.param)
       
     }
     c0<-matrix(0,nn,n)
     for ( i in 1: nn) 
      for ( j in 1:n) {  
		        sdist<-sqrt((floc[i,1]-loc[j,1])^2+(floc[i,2]-loc[j,2])^2)           
			tdist<-abs(ftimes[i]-times[j])	
         		c0[i,j]<-cov.theo(sdist,tdist,cov.param)
     }


     return(list(covz=covz,c0=c0))	     

}



