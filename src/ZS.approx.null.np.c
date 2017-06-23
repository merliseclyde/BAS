#include <R.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "sampling.h"

void LogBF_ZS_null_vect(double *r2curr, int *n, int *dim, int *nmodels, double *logmarg) {

  double LogBF_ZS_null(double r2curr, int n, int d);
  double ZS_logmarg(double r2curr, int n, int d, double rscale);
  double rscale = 1.0;

  int i;
     for (i=0; i < *nmodels; i++) {
//        logmarg[i] = LogBF_ZS_null(r2curr[i], *n, dim[i]);
          logmarg[i] = ZS_logmarg(r2curr[i], *n, dim[i], rscale);
     }
}

typedef struct C_int_struct
{
 void (*f)(double *x, int n, SEXP theta) ;    /* function */
 SEXP theta;  /*other args to f */
} C_int_struct, *IntStruct;

static void Cintfn(double *x, int n, void *ex)
{
  int i;

  IntStruct IS = (IntStruct) ex;

  IS->f(x, n, IS->theta);

    for( i=0; i<n; i++) {
    if(!R_FINITE(x[i]))
      error("non-finite function value");
  }
  return;
}

// use R's integrate code from QUADPACK to obtain marginal likelihood
double ZS_logmarg(double R2, int n, int d, double rscale) {

  double bound=0.0, epsabs, epsrel, result, abserr, *work, *ex;
  int inf = 1L, neval, ier, limit=200, lenw, last, *iwork;
  SEXP Rtheta;
  C_int_struct is;

  if (d <= 1 || n - d <= 1) return(0.0);

  epsabs = R_pow(DOUBLE_EPS, 0.25);
  epsrel = epsabs;
  lenw = 4 * limit;
  iwork = (int *) R_alloc((size_t) limit, sizeof(int));
  work = (double *) R_alloc((size_t) lenw, sizeof(double));

  PROTECT(Rtheta = allocVector(REALSXP, 4));
  REAL(Rtheta)[0] = R2;
  REAL(Rtheta)[1] = (double) n;
  REAL(Rtheta)[2] = (double) d;
  REAL(Rtheta)[3] = (double) rscale;

  ex = REAL(Rtheta);

  is.f = ZS_density;
  is.theta = Rtheta;

  Rdqagi(Cintfn, (void*)&is, &bound,&inf,&epsabs,&epsrel,&result,
         &abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);
//  Rprintf("ZS return: logBF %lf R2=%lf n= %lf d=%lf r=%lf \n", log(result), ex[0], ex[1], ex[2], ex[3]);

  UNPROTECT(1);
  return(log(result));
}

double ZS_shrinkage(double R2, int n, int d, double rscale) {

  double bound=0.0, epsabs, epsrel, result, abserr, *work, *ex;
  int inf = 1L, neval, ier, limit=200, lenw, last, *iwork;
  SEXP Rtheta;
  C_int_struct is;

  if (d <= 1) return(1.0);

  epsabs = R_pow(DOUBLE_EPS, 0.25);
  epsrel = epsabs;
  lenw = 4 * limit;
  iwork = (int *) R_alloc((size_t) limit, sizeof(int));
  work = (double *) R_alloc((size_t) lenw, sizeof(double));

  PROTECT(Rtheta = allocVector(REALSXP, 4));
  REAL(Rtheta)[0] = R2;
  REAL(Rtheta)[1] = (double) n;
  REAL(Rtheta)[2] = (double) d;
  REAL(Rtheta)[3] = (double) rscale;

  ex = REAL(Rtheta);

  is.f = ZS_density_shrinkage;
  is.theta = Rtheta;

  Rdqagi(Cintfn, (void*)&is, &bound,&inf,&epsabs,&epsrel,&result,
         &abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);


  result = exp(log(result) - ZS_logmarg(R2, n, d, rscale));
//  Rprintf("Shrinkage return: logBF %lf R2=%lf n= %lf d=%lf r=%lf \n", result, ex[0], ex[1], ex[2], ex[3]);
  UNPROTECT(1);
  return(result);
}

void ZS_density(double *x, int n, SEXP Rex) {
// d is p + 1  and includes the intercept
// prior for  1/g ~ gamma(1/2, rscale*n/2)
   double g, R2, rscale, d, nobs;
   int i;

   PROTECT(Rex);
   // SEXP Rex = PROTECT(duplicate(Rtheta));

   R2 =  REAL(Rex)[0];
   nobs = REAL(Rex)[1];
   d = REAL(Rex)[2];
   rscale = REAL(Rex)[3];
   for (i=0; i < n; i++) {
    g = x[i];
    x[i]  = .5*(log(1.0 + g)* (nobs - d) - (nobs-1.0)*log(1.0 + (1.0 - R2)*g));
    x[i] += .5*(log(.5*nobs*rscale) -3.0*log(g) - rscale* nobs/g) - lgamma(.5);

//    Rprintf("integrate: g= %lf BF= %lf R2=%lf nobs= %lf d=%lf r=%lf \n",
//            g, exp(x[i]), R2, nobs, d, rscale);

    x[i] = exp(x[i]);
  }
 UNPROTECT(1);
}

void ZS_density_shrinkage(double *x, int n, SEXP Rex) {
  // d is p + 1  and includes the intercept
  // prior for  1/g ~ gamma(1/2, rscale*n/2)
  double g, R2, rscale, d, nobs;
  int i;

  PROTECT(Rex);

  R2 =  REAL(Rex)[0];
  nobs = REAL(Rex)[1];
  d = REAL(Rex)[2];
  rscale = REAL(Rex)[3];
  for (i=0; i < n; i++) {
    g = x[i];
    x[i]  = .5*(log(1.0 + g)* (nobs - d) - (nobs-1.0)*log(1.0 + (1.0 - R2)*g));
    x[i] += .5*(log(.5*nobs*rscale) -3.0*log(g) - rscale* nobs/g) - lgamma(.5);
    x[i] = exp(x[i])*g/(1.0+g);
  }
  UNPROTECT(1);
}



double LogBF_ZS_null(double R2, int n, int d) {

/* this computes a Laplace approximation to the log of the Bayes factor
   with respect to the null model (intercept only), log(m_k)-log(m_0)

   R2 = 1-SSE/SST; the coefficient of determination
   n  = sample size;
   k  = number of covariates of the current model

   The prior under consideration is Zellner & Siow (1980)
*/

/* this version: JAN 14 2003 */

  void posroot(double a, double b, double c, double *root, double *status);
  double lik_null(double g, double R2,int n,int k);
  double info_null(double g, double R2, int n, int k);

  double status,root;
  double a,b,c,aux;
  int k;

  k = d - 1;
  aux=-(1-R2)*((double)k+3.);
  a=(double)n-4.-(double)k-2.*(1-R2);
  b=(double)n*(2.-R2)-3.;
  c=(double)n;

  a=a/aux;
  b=b/aux;
  c=c/aux;


  posroot(a,b,c,&root,&status);
    if (k == 0 || n <= k+1 || R2 >= 1.0) { return(0.0); }
    else {
      if(status!=1.){
	if(status==0.) {
	  Rprintf("No positive roots R2=%lf n=%d k=%d\n",R2,n,k);}
	else {
	  Rprintf("\n More than one positive root  R2=%lf n=%d k=%d\n",R2,n,k);}
      }
      //    else{return(lik_null(root,R2,n,k)+(log(4.*asin(1.)) - log(info_null(root,R2,n,k)))/2.);}
      else{return(lik_null(root,R2,n,k)+ log(sqrt(2.*PI)) - .5*log(info_null(root,R2,n,k)));}
    }
  return(NA_REAL);
}

double lik_null(double g, double R2, int n, int k){
/* this computes log(likelihood  x prior), where the likelihood is marginal
   on the intercept, regression coefficients and the variance
*/
  double aux;
  if (R2 >= 1.0) R2 = 1.0;
  aux=((double)n-1.-(double)k)*log(1.+g)-((double)n-1.)*log(1.+(1.-R2)*g)-
    3.*log(g)-((double)n)/g;
  aux=aux/2.;
  aux=aux+log((double)n/2.)/2.-lgammafn(.5);
  return(aux);
}


double info_null(double g,double R2,int n,int k){
  double aux;

  aux= -((double)n-1.-(double)k)/R_pow_di(1.+g,2);
  aux=aux+((double)n-1.)*R_pow_di(1.-R2,2)/R_pow_di(1.+(1.-R2)*g,2)+3/R_pow_di(g,2);
  aux=aux-2.*(double)n/R_pow_di(g,3);
  aux=-aux/2.;
  return(aux);
}

void posroot(double a, double b, double c, double *root, double *status)
{ /* this computes the real roots of a cubic polynomial; in the end, if
     status==1, root stores the nonegative root; if status is not one,
     then status is the total number of nonegative roots and root is
     useless
  */

  int i;
  double Q,R,disc,Q3,A,B,aux,x[3];

  *root = 0.;
  *status=0.;

  Q=(R_pow_di(a,2)-3.*b)/9.;
  R=(2*R_pow_di(a,3)-9*a*b+27.*c)/54.;
  Q3=R_pow_di(Q,3);

  disc=R_pow_di(R,2)-Q3;

  if(disc>=0.){
    if(R>=0) A=-cbrt(R+sqrt(disc));
    else A=-cbrt(R-sqrt(disc));

    if(A==0.) B=0.;
    else B=Q/(A);
    *root=(A+B)-a/3.;
    if(*root>=0) *status=1.;
  }
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
}
