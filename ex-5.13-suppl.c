#include <math.h>
#include <R.h>
#include <Rmath.h>
#define INF 10.0e32


void matneigh(double *x, double *y, int *n, int *m, int *neighborhood_type)
{
/* 
   n number of rows
   m number of columns
   neighbor_type
   case 0 :  4 nearest neighborhood
   case 1 :  8 nearest neighborhood
   case 2 : 12 nearest neighborhood
   
 */  
	int  edge, i, j, k, l,nm;
	edge=1;
        if (*neighborhood_type > 1) 
		    edge=2;

	k=0;
	nm=(*n-2*(edge))*(*m-2*(edge));

	for (i=(0+(edge)); i <(*n-(edge)); i++) {			
	    for (j=(0+(edge)); j <(*m-(edge)); j++) {		

	        l=0;
	    	if(*neighborhood_type == 0) {
			    y[k+l*(nm)]=(x[i-1+j*(*n)]+x[i+1+j*(*n)]);
			    l=l+1;
			    y[k+l*(nm)]=(x[i+(j-1)*(*n)]+x[i+(j+1)*(*n)]);		
			}
		if (*neighborhood_type == 1) {
			    y[k+l*(nm)]=(x[i-1+j*(*n)]+x[i+1+j*(*n)]);
			    l=l+1;
			    y[k+l*(nm)]=(x[i+(j-1)*(*n)]+x[i+(j+1)*(*n)]);
			    l=l+1;
			    
			    y[k+l*(nm)]=(x[i+1+(j+1)*(*n)]+x[i-1+(j-1)*(*n)]);
			    l=l+1;
			    y[k+l*(nm)]=(x[i+1+(j-1)*(*n)]+x[i-1+(j+1)*(*n)]);
			}
		if (*neighborhood_type == 2) {
			    y[k+l*(nm)]=(x[i-1+j*(*n)]+x[i+1+j*(*n)]);
			    l=l+1;
			    y[k+l*(nm)]=(x[i+(j-1)*(*n)]+x[i+(j+1)*(*n)]);
			    l=l+1;
			    
			    y[k+l*(nm)]=(x[i+1+(j+1)*(*n)]+x[i-1+(j-1)*(*n)]);
			    l=l+1;
			    y[k+l*(nm)]=(x[i+1+(j-1)*(*n)]+x[i-1+(j+1)*(*n)]);
			    l=l+1;
			    
			    y[k+l*(nm)]=(x[i-2+j*(*n)]+x[i+2+j*(*n)]);
			    l=l+1;
			    y[k+l*(nm)]=(x[i+(j-2)*(*n)]+x[i+(j+2)*(*n)]);
		}
		k=k+1;	
	    }
	    
	}	
	
/*	PutRNGstate();	*/
}


void latcov(double *x, double *xcov, int *n, int *m, int *dimilag, int *dimjlag)
{
/* 
   n number of rows
   m number of columns
   neighbor_type
   case 0 :  4 nearest neighborhood
   case 1 :  8 nearest neighborhood
   case 2 : 12 nearest neighborhood
   
 */  
	int i, j, k, l;
	double mu,temp;
	mu=0;
        for (i=0; i <(*n); i++) {
	    for (j=0; j <(*m); j++) {
		mu+=x[i+j*(*n)];	
	    }
	}	
	mu=mu/(*n*(*m));

	   for (k=0; k<(*dimilag); k++) {
	    for (l=0; l <(*dimjlag); l++) {
	        xcov[k+l*(*dimilag)]=0;
		    for (i=(*n-1); i>=(k); i--) {
			for (j=0; j <(*m-l); j++) {
			    xcov[k+l*(*dimilag)]+=(x[i+j*(*n)]-mu)*(x[i-k+(j+l)*(*n)]-mu);
			}    
		    }
		    xcov[k+l*(*dimilag)]=xcov[k+l*(*dimilag)]/((*n)*(*m));       		
	   }
	}	   
	
}

void one_gibbs(double *x, int *n, int *m, double *beta, int *neighborhood_type,int  edge)
{
/* 
   n number of rows
   m number of columns
   neighbor_type
   case 0 :  4 nearest neighborhood
   case 1 :  8 nearest neighborhood
   case 2 : 12 nearest neighborhood
   
 */  
	int  i, j, k;
	double temp, u;
	for (i=(0+edge); i <(*n-edge); i++) {		
	    for (j=(0+edge); j <(*m-edge); j++) {		
	    
	    	if(*neighborhood_type == 0) {
			    temp=beta[0]+beta[1]*(x[i-1+j*(*n)]+x[i+1+j*(*n)]);
			    temp+=beta[2]*(x[i+(j-1)*(*n)]+x[i+(j+1)*(*n)]);
			}
		else if (*neighborhood_type == 1) {
			    temp=beta[0]+beta[1]*(x[i-1+j*(*n)]+x[i+1+j*(*n)]);
			    temp+=beta[2]*(x[i+(j-1)*(*n)]+x[i+(j+1)*(*n)]);
			    temp+=beta[3]*(x[i+1+(j+1)*(*n)]+x[i-1+(j-1)*(*n)]);
			    temp+=beta[4]*(x[i+1+(j-1)*(*n)]+x[i-1+(j+1)*(*n)]);
			}
		else {
			    temp=beta[0]+beta[1]*(x[i-1+j*(*n)]+x[i+1+j*(*n)]);
			    temp+=beta[2]*(x[i+(j-1)*(*n)]+x[i+(j+1)*(*n)]);
			    temp+=beta[3]*(x[i+1+(j+1)*(*n)]+x[i-1+(j-1)*(*n)]);
			    temp+=beta[4]*(x[i+1+(j-1)*(*n)]+x[i-1+(j+1)*(*n)]);
			    temp+=beta[5]*(x[i-2+j*(*n)]+x[i+2+j*(*n)]);
			    temp+=beta[6]*(x[i+(j-2)*(*n)]+x[i+(j+2)*(*n)]);
		}


	    
	    
	    
    	    u=unif_rand();	    
	    temp=1/(1+exp(temp));
	
	    if (u<= temp) {
		x[i+j*(*n)]=0.0;
		}
	    else {
		x[i+j*(*n)]=1.0;
		}		

	    }
	}	
					
	
/*	PutRNGstate();	*/
}

void gibbs(double *x, int *n, int *m, double *beta, int *neighborhood_type, int *burn_in, int *iter)
{
/*
burn_in number of the iterations for the burning
iter    number of the iterations for the MCMC algorithm
*/
    	int i,k, edge;
	GetRNGstate(); 
	Rprintf("Burn-in \n");
        edge=1;
        if (*neighborhood_type > 1) 
		    edge=2;
	
	for (i=0; i< *burn_in; i++) {
/*	    Rprintf("%d \n",i+1); */
	    one_gibbs(x, n, m, beta, neighborhood_type, edge);    
	}
	Rprintf("MCMC \n");
    	for (i=0; i< *iter; i++) {

	    one_gibbs(x, n, m, beta, neighborhood_type, edge);	  
	}
	PutRNGstate();	
}


double potential(int i1, int j1, double *x, int *n, double *beta, int *neighborhood_type)
{
	double temp;
	int i,j;
	i=i1;
	j=j1;


    		if(*neighborhood_type == 0) {
			    temp=beta[0]+beta[1]*(x[i-1+j*(*n)]+x[i+1+j*(*n)]);
			    temp+=beta[2]*(x[i+(j-1)*(*n)]+x[i+(j+1)*(*n)]);
			}
		else if (*neighborhood_type == 1) {
			    temp=beta[0]+beta[1]*(x[i-1+j*(*n)]+x[i+1+j*(*n)]);
			    temp+=beta[2]*(x[i+(j-1)*(*n)]+x[i+(j+1)*(*n)]);
			    temp+=beta[3]*(x[i+1+(j+1)*(*n)]+x[i-1+(j-1)*(*n)]);
			    temp+=beta[4]*(x[i+1+(j-1)*(*n)]+x[i-1+(j+1)*(*n)]);
			}
		else {
			    temp=beta[0]+beta[1]*(x[i-1+j*(*n)]+x[i+1+j*(*n)]);
			    temp+=beta[2]*(x[i+(j-1)*(*n)]+x[i+(j+1)*(*n)]);
			    temp+=beta[3]*(x[i+1+(j+1)*(*n)]+x[i-1+(j-1)*(*n)]);
			    temp+=beta[4]*(x[i+1+(j-1)*(*n)]+x[i-1+(j+1)*(*n)]);
			    temp+=beta[5]*(x[i-2+j*(*n)]+x[i+2+j*(*n)]);
			    temp+=beta[6]*(x[i+(j-2)*(*n)]+x[i+(j+2)*(*n)]);
		}

		
	return temp;
		
}

void spinexchange(double *x, int *n, int *m, double *beta, int *neighborhood_type,  int *iter)
{
/*
burn_in number of the iterations for the burning
iter    number of the iterations for the MCMC algorithm
*/
    	int k, edge,i,i1,j,j1,m1,n1;
	double a,a1,temp,u;
        edge=1;
        if (*neighborhood_type > 1) 
		    edge=2;
	GetRNGstate(); 
/*	Rprintf("Burn-in \n"); */
	for (k=0; k< *iter; k++) {
	    temp=(*n-1)*unif_rand();
	    i=floor(temp)+edge;

	    temp=(*m-1)*unif_rand();
	    j=floor(temp)+edge; 
	    temp=(*n-1)*unif_rand();
	    i1=floor(temp)+edge; 
	    temp=(*m-1)*unif_rand();
	    j1=floor(temp)+edge; 
	    
	    *n=*n+2*edge;
            *m=*m+2*edge;


	    a=potential(i, j, x,  m, beta, neighborhood_type);
	    a1=potential(i1, j1, x, m, beta, neighborhood_type);
            temp=(x[i+(j)*(*n)]-x[i1+(j1)*(*n)])*a+(x[i1+(j1)*(*n)]-x[i+(j)*(*n)])*a1;
	   
/* exchange */	    
	    if (temp > 0) {
	           u=x[i+(j)*(*n)];
		   x[i+(j)*(*n)]=x[i1+(j1)*(*n)];
		   x[i1+(j1)*(*n)]=u;
		    	}
	    else {
	        u=unif_rand();
	        if ( log(u) < temp ) {
		   u=x[i+(j)*(*n)];
		   x[i+(j)*(*n)]=x[i1+(j1)*(*n)];
		   x[i1+(j1)*(*n)]=u;
	        }
	    }	
/*	    Rprintf("%d \n",i+1); */

	    *n=*n-2*edge;
            *m=*m-2*edge;
	
	}

	PutRNGstate();	
}
