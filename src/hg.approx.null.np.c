#include "sampling.h"
#define lgamma lgammafn
double LogBF_Hg_null(double r2curr, int n, int d, double alpha, int gpower);
void posroot(double a, double b, double c, double *root, double *status);
double lik_null_HG(double g, double R2,int n,int k, double alpha, int gpower);
double info_null_HG(double g, double R2, int n, int k, double alpha);

double lik_null_HG(double g, double R2, int n, int k, double alpha, int gpower){
/* this computes log(likelihood  x prior), where the likelihood is marginal
   on the intercept, regression coefficients and the variance
*/
  double aux;

  aux=((double)n-1.-(double)k)*log(1.+g)-((double)n-1.)*log(1.+(1.-R2)*g)+ 2.*(double)gpower*log(g)-alpha*log(1.+g/(double)n);
  aux=aux/2.;
  aux=aux - log((double)n) +  log(alpha/2.-1.); 

  return(aux);
}


double info_null_HG(double g,double R2,int n,int k, double alpha)
{/* This computs the second derivative of LogLik(tau) which is
    equal to
 (1/2)*(- alpha ng/(n+g)^2-(n-1) eg/(1+eg)^2+(n-1-p)g/(1+g)^2)
where g=e^{\tau}
*/

  double aux;

  aux= ((double)n-1.-(double)k)*g/R_pow_di(1.+g,2);
  aux=aux-((double)n-1.)*(1.-R2)*g/R_pow_di(1.+(1.-R2)*g,2);
  aux=aux-alpha*(double)n*g/R_pow_di((double)n+g,2);
  aux=aux/2.;
  return(aux);
}



double LogBF_Hg_null(double R2, int n, int d, double alpha, int gpower){

/* this computes a Laplace approximation to the log of the Bayes factor
   with respect to the null model (intercept only), log(m_k)-log(m_0)

   R2 = 1-SSE/SST; the coefficient of determination
   e = 1 - R2; 
   n  = sample size;
   k  = number of covariates of the current model (exclusing the intercept)

   The prior under consideration is Hyper-g with prior (1+g/n)^{-alpha/2}
   
   The cubic equation for loglik'(g)=0 is
   -e (alpha - 2 gpower + k) g^3 
   - {[e (k - 2 gpower) - (1 - e)] n + (1 - e) + k + (1 + e) (alpha - 2 gpower)} g^2 
   + [(1 - e) (n - 1) n - k n - alpha + 2 gpower (1 + (1 + e) n)] g + 2 gpower n
*/

/* this version: April 11 2005 */


  double status,root, logmarg;
  double a,b,c,e,aux;
  int k;

  logmarg = NA_REAL;
  k = d - 1;
  e = 1.-R2;
  aux=-e*(alpha - (double)gpower*2. + (double)k);
  a= - ((e*((double)k - 2.*(double)gpower) - (1. - e))*(double)n + (1. - e) + (double)k + (1. + e)*(alpha - 2.*(double)gpower));
  b= ((1. - e)*((double)n - 1.)*(double)n - (double)k*(double)n - alpha + 2.*(double)gpower*(1. + (1. + e)*(double)n));
  c=(double)n*2.*(double)gpower;

  a=a/aux;
  b=b/aux;
  c=c/aux;
  posroot(a,b,c,&root,&status);
  if (k == 0 || n <= d) { logmarg = 0.0; }
  else {
    if(status!=1.){
      if(status==0.) Rprintf("\n No positive roots\n");
      else Rprintf("\n More than one positive root; this should not happen\n");
    }
    else{
      logmarg = lik_null_HG(root,R2,n,k, alpha, gpower)+
	(log(4.*asin(1.))-
	 log(-info_null_HG(root,R2,n,k, alpha)))/2.;}
  }
  return(logmarg);

}


void LogBF_Hg_null_vect(double *r2curr, int *n, int *dim, int *nmodels, double *logmarg, double *alpha, int *gpower) {


  int i;
     for (i=0; i < *nmodels; i++) {
        logmarg[i] = LogBF_Hg_null(r2curr[i], *n, dim[i], *alpha, *gpower);
     }
}
