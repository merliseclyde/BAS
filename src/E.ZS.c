#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#define lgamma lgammafn
double E_ZS_approx_null(double R2, int n, int k){

/* this computes a Laplace approximation to the posterior expectation
   of g/(1+g) under model M_k with $k$ covariates; there is a common
   intercept

   E[g/(1+g) | y, M_k]

   R2 = 1-SSE/SST; the coefficient of determination
   n  = sample size;
   k  = number of covariates of the current model

   The prior under consideration is Zellner & Siow (1980)
*/

/* this version: JAN 19 2003 */

  void posroot(double a, double b, double c, double *root, double *status);
  double h1(double g, double eps, int n, int k);
  double h2(double g, double eps, int n, int k);
  double infoh1(double g, double eps, int n, int k);
  double infoh2(double g, double eps, int n, int k);

  double status,root1,root2;
  double a,b,c,aux,eps,denom,numer, s;

  eps=1.0-R2;

  aux=-eps*((double)k+3.);
  a=(double)n-(double)k-4.;
  b=(double)n*(1.+eps)-1.;
  c=(double)n;

  a=a/aux;
  b=b/aux;
  c=c/aux;
  numer = 1.0;
  denom = 1.0;

  if (k == 0 || n <= k+1 || R2 >= 1.0) { return(0.0); }
  else {
    posroot(a,b,c,&root1,&status);
  
   if(status!=1.){
     if(status==0.) error("\n No positive roots for the numerator  R2=%lf n=%d k=%d\n\n");
     else error("\n More than one positive root for the numerator\n");
   }
   else{
     numer=h1(root1,eps,n,k)-log(-infoh1(root1,eps,n,k))/2.;
   }

   aux=-eps*((double)k+3.);
   a=(double)n-(double)k-2.*eps-4.;
   b=(double)n*(1.+eps)-3.;
   c=(double)n;

   a=a/aux;
   b=b/aux;
   c=c/aux;


   posroot(a,b,c,&root2,&status);

   if(status!=1.){
     if(status==0.) error("\n No positive roots for the denominator  R2=%lf n=%d k=%d\n\n");
     else error("\n More than one positive root for the denominator\n");
   }
   else{
     denom=h2(root2,eps,n,k)-log(-infoh2(root2,eps,n,k))/2.;
   }
   s = exp(numer-denom);
   return(s);
  }
}

double h1(double g, double eps, int n, int k){
  double aux;

  aux=((double)n-3.-(double)k)*log(1.+g)-((double)n-1.)*log(1.+eps*g)-
    log(g)-((double)n)/g;
  aux=aux/2.;
  return(aux);
}

double h2(double g, double eps, int n, int k){
  double aux;

  aux=((double)n-1.-(double)k)*log(1.+g)-((double)n-1.)*log(1.+eps*g)-
    3*log(g)-((double)n)/g;
  aux=aux/2.;
  return(aux);
}

double infoh1(double g, double eps, int n, int k){
  double aux;
  aux= -((double)n-3.-(double)k)/R_pow_di(1.+g,2);
  aux=aux+((double)n-1.)*R_pow_di(eps,2)/R_pow_di(1.+eps*g,2)+1./R_pow_di(g,2);
  aux=aux-2.*(double)n/R_pow_di(g,3);
  aux=aux/2.;
  return(aux);
}

double infoh2(double g, double eps, int n, int k){
  double aux;
  aux= -((double)n-1.-(double)k)/R_pow_di(1.+g,2);
  aux=aux+((double)n-1.)*R_pow_di(eps,2)/R_pow_di(1.+eps*g,2)+3./R_pow_di(g,2);
  aux=aux-2.*(double)n/R_pow_di(g,3);
  aux=aux/2.;
  return(aux);
}

