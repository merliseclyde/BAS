#include <R.h>
#include <math.h>
#include <Rmath.h>

/*  vectorized version of Log BFZS full that can be called from R directly */

 void posroot_full(double a, double b, double c, double *root, double *status);
 double lik_full(double g, double eps ,int n, int p, int k);
 double info_full(double g, double eps, int n, int p, int k);

 // legacy code not in use currently
 // start no cov
void LogBF_ZS_full_vect(double *R2full, double *r2curr, int *n, int *ptotal, int *dim, int *nmodels, double *logmarg) {
  double LogBF_ZS_full(double r2full, double r2curr, int n, int ptotal, int d);

  int i;
     for (i=0; i < *nmodels; i++) {
        logmarg[i] = LogBF_ZS_full(*R2full , r2curr[i], *n, *ptotal, dim[i]);
     }
}
// end no cov

double LogBF_ZS_full(double r2full, double r2curr, int n, int ptotal, int d){
/*
   Zellner-Siow prior for the linear model; here we compute the Bayes
   factor to the full model using a Laplace approximation.
  
   The current model has $k$ covariates = d - 1  of a total of $p$ (ptotal - 1)    possible and  
   there is a common intercept in all models being compared. The sample 
   size is $n$.
   
   Originally we are considering
  
     Y= 1 alpha + X1 beta1 + X2 beta2
   
   where 1'X1=1'X2=0, rank(X1)=k, rank(X2)=p-k, and we are interested in
   testing H0: beta2=0.
  
   We then reparametrize to
  
     Y= 1 alpha + X1 eta + V beta2
  
   where eta = beta1 + (X1'X1)^-1 X2 beta2
           V = [I - X1(X1'X1)^-1 X1'] X2
   
                 y'(I-P*-P_X1-P_V)y       y'(I-P*-P_X1-P_X2)y
   R^2Full =  1- ------------------- = 1- -------------------
                      y'(I-P*)y              y'(I-P*)y
  
                    y'(I-P*-P_X1)y
   R^2current =  1- --------------
                      y'(I-P*)y
  
     P* = 1 (1'1)^(-1) 1'= 1 1'/n, where 1 is the nx1 vector of ones
     P_X1  = X1 (X1'X1)^(-1) X1'
     P_V  = V (V'V)^(-1) V'
     I  = nxn identity matrix

   THIS VERSION: JAN 14 2003
*/


  double status,root, ZSlogmarg=0.0;
  double eps,a,b,c,aux;
  int k, p;

  p = ptotal - 1;
  k = d - 1;
  eps=(1.-r2full)/(1.-r2curr);
  aux=-eps*((double)p-(double)k+3.);
  a=(double)n-(double)p+eps*((double)k-2.)-4.;
  b=(double)n*(eps+1.)-3.;
  c=(double)n;

  a=a/aux;
  b=b/aux;
  c=c/aux;

  posroot_full(a,b,c,&root,&status);

  // start no cov
  if(status!=1.){
    if(status==0.) Rprintf("\n No positive roots\n");
    else Rprintf("\n More than one positive root\n");
  } // end  no. cov
  else{
    if (k != p) {
      ZSlogmarg = -lik_full(root,eps,n,p,k)+(-log(4.*asin(1.))+ log(-info_full(root,eps,n,p,k)))/2.;
    }
    else   ZSlogmarg = 0.0;
  }
  return(ZSlogmarg);
}

double lik_full(double g, double eps, int n, int p, int k){
  double aux;
  
  aux=((double)n-((double)p+1.))*log(1.+g);
  aux=aux-((double)n-((double)k+1.))*log(1.+eps*g)-
    3.*log(g)-((double)n)/g;
  aux=aux/2.;
  aux=aux+log((double)n/2.)/2.-lgammafn(.5);
  return(aux);
}


double info_full(double g, double eps, int n, int p, int k){
  double aux;

  aux=-((double)n-((double)p+1.))/R_pow_di(1.+g,2);
  aux=aux+((double)n-((double)k+1.))*R_pow_di(eps,2)/R_pow_di(1.+eps*g,2);
  aux=aux+3./R_pow_di(g,2)-2*(double)n/R_pow_di(g,3);
  aux=aux/2.;
  return(aux);

  return(aux);
}

void posroot_full(double a, double b, double c, double *root, double *status)
{
  int i;
  double Q,R,disc,Q3,A,B,aux,x[3];

  *status=0.;

  Q=(R_pow_di(a,2)-3.*b)/9.;
  R=(2*R_pow_di(a,3)-9*a*b+27.*c)/54.;
  Q3=R_pow_di(Q,3);

  disc=R_pow_di(R,2)-Q3;

  if(disc>=0.){
    if(R>=0) A=-cbrt(R+sqrt(disc));
    else A=-cbrt(R-sqrt(disc));

    if(A==0.)B=0.;
    else B=Q/A;
    *root=(A+B)-a/3.;
    if(*root>=0) *status=1.;
  }
  // start no cov
  // negative root should not occur
  else{
    A=acos(R/sqrt(Q3));
    aux= 2. * sqrt(Q);
    x[0]=-aux * cos(A/3.);
    x[1]=-aux * cos((A+4.*asin(1.))/3.);
    x[2]=-aux * cos((A-4.*asin(1.))/3.);
    aux=a/3.;
    for(i=0;i<3;i++) x[i]=x[i]-aux;
    for(i=0;i<3;i++){
      if (x[i]>=0.){
      	*status=*status+1.;
	      *root=x[i];
      }
    }
  }
  // end no cov
}

