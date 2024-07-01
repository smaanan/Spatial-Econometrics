#include <stdio.h>
#include <stdlib.h> /* for malloc & free */
#include "memory.h"
#include <math.h> 
#include <R.h>
#include <Rmath.h>

#define Integer int
#define Real double

void empvar(Integer *n, Real *xc, Real *yc, Real *tc, Real *data, 
	   Integer *nbins, Integer *ntbins, Real *lims, Real *tlims, Integer *modulus, 
	    Real *maxdist, Real *tmaxdist, Integer *cbin, Real *vbin,
	   Integer *sdcalc, Real *sdbin)
/*
n number of observations
nbins number of spatial bins
ntbins number of temporal bins
xc x coordinates
yc y coordinates
tc t coordinates
lims limits of the interval for the spatial components
tlims limits of the interval for the temporal components
maxdist maximum spatial distance
tmaxdist maximum time distance
vbins mean of observations in a bin (semi-variogram estimate)  vector of nbins x ntbins
      the order is v(ds_1,dt_1),v(ds_2,dt_1),...v(ds_nbins,dt_1),...,v(ds_1,dt_ntbins,),....v(ds_nbins,dt_ntbins)
cbins count the number of points in a bin
*/
{
  
  Integer i, j, ind=0, tind=0;

  Real v=0.0;
  Real dist=0.0, dx=0.0, dy=0.0, dt=0.0, tdist=0.0;

  for (j=0; j < *n; j++)
    { 
      for (i=j+1; i<*n; i++) 
	{
	  dx = xc[i] - xc[j];
	  dy = yc[i] - yc[j];
          dt = tc[i] - tc[j];
	  dist = pythag(dx, dy);
	  tdist = fabs(dt);
	  
	  if((dist <= *maxdist) && (tdist <= *tmaxdist) )
	    {
	      v = data[i] - data[j];
	      if (*modulus) v = sqrt(sqrt(v*v));
	      else v = (v*v)/2.0;
	      ind = 0;
	      tind = 0;

	      while (dist > lims[ind] && ind <= *nbins ) ind++ ;

	      while (tdist > tlims[tind] && tind <= *ntbins ) tind++ ;	
	      if ((dist <= lims[ind]) && (tdist <= tlims[tind]))
		{
		  vbin[(ind-1+(tind-1)*(*nbins))]+= v; 
		  cbin[(ind-1+(tind-1)*(*nbins))]++;
		  if(*sdcalc) sdbin[(ind+(tind-1)*(*nbins)-1)] += v*v;
		}

	    }
	}
    }
  
  for (j=0; j < (*nbins*(*ntbins)); j++) 
    {
      if (cbin[j]){
	if(*sdcalc)
	  { 
	    sdbin[j] = sqrt((sdbin[j] - ((vbin[j] * vbin[j])/cbin[j]))/(cbin[j] - 1));
	  }
	vbin[j] = vbin[j]/cbin[j];
	if (*modulus) {
	  vbin[j] = vbin[j] * vbin[j];
	  vbin[j] = (vbin[j] * vbin[j])/(0.914 + (0.988/cbin[j]));
	}
      }
    }
  
}

