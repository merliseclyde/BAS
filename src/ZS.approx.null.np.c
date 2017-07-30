#include <R.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "bas.h"

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
//    if(!R_FINITE(x[i]))
//      Rprintf("warning: non-finite function value for integral\n");
  }
  return;
}

// use R's integrate code from QUADPACK to obtain marginal likelihood
double ZS_logmarg(double R2, int n, int d, double rscale) {

  double mode = 0.0, status=1.0, root, bound=DBL_EPSILON, epsabs, epsrel, result, abserr, *work, *ex;
  int inf = 1L, neval, ier, limit=200, lenw, last, *iwork;
  SEXP Rtheta;
  C_int_struct is;



  if (d <= 1 || n - d <= 1) return(0.0);

  epsabs = R_pow(DOUBLE_EPS, 0.25);
  epsrel = epsabs;
  lenw = 4 * limit;
  iwork = (int *) R_alloc((size_t) limit, sizeof(int));
  work = (double *) R_alloc((size_t) lenw, sizeof(double));

  mode = find_mode_g_JZS(R2, n, d, &root, &status);  // log density at mode
 // Rprintf("mode %lf g-mode %lf status = %lf\n", mode, root, status);

  PROTECT(Rtheta = allocVector(REALSXP, 6));
  REAL(Rtheta)[0] = R2;
  REAL(Rtheta)[1] = (double) n;
  REAL(Rtheta)[2] = (double) d;
  REAL(Rtheta)[3] = (double) rscale;
  REAL(Rtheta)[4] = (double) mode;
  REAL(Rtheta)[5] = (double) root;

  ex = REAL(Rtheta);

  is.f = ZS_density;
  is.theta = Rtheta;

  inf = 2L;  // use u = log(g) os integration range i -Inf to Inf

  Rdqagi(Cintfn, (void*)&is, &bound,&inf,&epsabs,&epsrel,&result,
         &abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);

//  if (!R_FINITE(result)) {
//    Rprintf("ZS return: logBF %lf R2=%lf n= %lf d=%lf r=%lf \n", log(result), ex[0], ex[1], ex[2], ex[3]);
//  }

  UNPROTECT(1);
  return(log(result) + mode);
}

double ZS_shrinkage(double R2, int n, int d, double rscale) {

  double mode, logmarg, bound=0.0, epsabs, epsrel, result, abserr, *work, *ex;
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
  REAL(Rtheta)[4] = (double) mode;
  ex = REAL(Rtheta);

  is.f = ZS_density_shrinkage;
  is.theta = Rtheta;

  Rdqagi(Cintfn, (void*)&is, &bound,&inf,&epsabs,&epsrel,&result,
         &abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);

  logmarg = ZS_logmarg(R2, n, d, rscale);
  if (!R_FINITE(result) | !R_FINITE(logmarg)) result = 1.0;
  else result = exp(log(result) - logmarg + mode);

//  Rprintf("Shrinkage return: logBF %lf R2=%lf n= %lf d=%lf r=%lf \n", result, ex[0], ex[1], ex[2], ex[3]);
  UNPROTECT(1);
  return(result);
}

// // return log(1 + exp(x)), preventing cancellation and overflow */
// From http://www.codeproject.com/KB/recipes/avoiding_overflow.aspx
// from  BayesFactor  bfutility
double log1pExp(double x)
{
  const double LOG_DBL_EPSILON = log(DBL_EPSILON);
  const double LOG_ONE_QUARTER = log(0.25);

  if (x > -LOG_DBL_EPSILON)
  {
    // log(exp(x) + 1) == x to machine precision
    return x;
  }
  else if (x > LOG_ONE_QUARTER)
  {
    return log( 1.0 + exp(x) );
  }
  else
  {
    // Prevent loss of precision that would result from adding small argument to 1.
    return log1p( exp(x) );
  }
}

double logExpXplusExpY( const double x, const double y )
{
  return x + log1pExp( y - x );
}

void ZS_density(double *x, int n, SEXP Rtheta) {
// d is p + 1  and includes the intercept
// prior for  1/g ~ gamma(1/2, rscale*n/2)
// x = log(g) for increased numerical stability
   double u, R2, rscale, d, nobs, mode, shift;
   int i;

   PROTECT(Rtheta);
   SEXP Rex = PROTECT(duplicate(Rtheta));

   R2 =  REAL(Rex)[0];

   nobs = REAL(Rex)[1];
   d = REAL(Rex)[2];
   rscale = REAL(Rex)[3];
   mode = REAL(Rex)[4];
   shift = REAL(Rex)[5];

   for (i=0; i < n; i++) {
    u = x[i] +  log(shift); // on input  log(g) = u; shift by mode
    x[i] = -mode;   //  subtract mode of f on log scale
    x[i] += .5*(log1pExp(u)*(nobs - d) -
               (nobs-1.0)*log1pExp(u + log(1.0 - R2)));
    x[i] += .5*log(.5*nobs*rscale) -lgamma(.5) -
            (1.0 + .5)*u - .5*rscale*nobs*exp(-u) ;
    x[i] += u;   // log jacobian
//  if (!R_FINITE(x[i]))  Rprintf("integrate: u= %lf logBF= %lf R2=%lf nobs= %lf d=%lf r=%lf mode=%lf shift%lf\n",
//            u, x[i], R2, nobs, d, rscale, mode, shift);

    x[i] = exp(x[i]);
  }
 UNPROTECT(2);
}

void ZS_density_g(double *x, int n, SEXP Rtheta) {
  // d is p + 1  and includes the intercept
  // prior for  1/g ~ gamma(1/2, rscale*n/2)
  double g, R2, rscale, d, nobs, mode;
  int i;

  PROTECT(Rtheta);
  SEXP Rex = PROTECT(duplicate(Rtheta));

  R2 =  REAL(Rex)[0];

  nobs = REAL(Rex)[1];
  d = REAL(Rex)[2];
  rscale = REAL(Rex)[3];
  mode = REAL(Rex)[4];

  for (i=0; i < n; i++) {
    g = x[i]; // on input
    x[i] = -mode;   // "divide" by mode or subtract on log scale
    x[i] += .5*(log(1.0 + g)*(nobs - d) - (nobs-1.0)*log(1.0 + (1.0 - R2)*g));
    x[i] += .5*(log(.5*nobs*rscale) -3.0*log(g) - rscale* nobs/g) - lgamma(.5);
    if (!R_FINITE(x[i]))  Rprintf("integrate: g= %lf BF= %lf R2=%lf nobs= %lf d=%lf r=%lf mode=%lf\n",
        g, exp(x[i]), R2, nobs, d, rscale, mode);

    x[i] = exp(x[i]);
  }
  UNPROTECT(2);
}

void ZS_density_shrinkage(double *x, int n, SEXP Rex) {
  // d is p + 1  and includes the intercept
  // prior for  1/g ~ gamma(1/2, rscale*n/2)
  double g, R2, rscale, d, nobs, mode;
  int i;

  PROTECT(Rex);

  R2     =  REAL(Rex)[0];
  nobs   = REAL(Rex)[1];
  d      = REAL(Rex)[2];
  rscale = REAL(Rex)[3];
  mode   = REAL(Rex)[4];

  for (i=0; i < n; i++) {
    g = x[i];  // on input
    x[i]  = -mode;
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

double find_mode_g_JZS(double R2, int n, int d, double *root, double *status) {
  int k;
  double aux, a, b, c, mode;
// d is dimension of model with intercept

void posroot(double a, double b, double c, double *root, double *status);
double lik_null(double g, double R2,int n,int k);

k = d - 1;
aux=-(1.0-R2)*((double)k+3.);
a=(double)n-4.-(double)k-2.*(1.0-R2);
b=(double)n*(2.-R2)-3.;
c=(double)n;

a=a/aux;
b=b/aux;
c=c/aux;

posroot(a,b,c,root,status);


if (*status != 1.0)  *root = n/20.0;  //should not occur but fix

mode = lik_null(*root,R2,n,k);  // this is log marginal at mode

return(mode);
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
