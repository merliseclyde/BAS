/*--------------------------------------------------------------------------
  CONHYPTWO1(a,b,c,x,y)
  Phi1 hypergeometric function of two variables.
  Gradshteyn and Ryzhik (1965), equation 9.261.1
  This program assumes:  0<a<c, 0<b, y<=1.
  x is a column vector of real numbers, other parameters scalar.
--------------------------------------------------------------------------*/
#include <R.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "bas.h"
#define FACCURACY 1e-12

/*--------------------------------------------------------------------------
  HYPERG2F1
  Function for the Hypergeometric2F1(a,b,c,z).
  Covers cases:
    {0<b<c, 0<=z<=1} or {0<b, abs(c-a-b)<1, z<0}.
--------------------------------------------------------------------------*/
double hyperg2F1(double a, double b, double c, double z) {
  int j=0;
  double w=1.0, F=1.0;
  if (a<0) {
    /* Apply Abramowitz and Stegun 15.3.3 */
    F=R_pow(1-z,c-a-b)*hyperg2F1(c-a,c-b,c,z);
  }
  else if (z<0) {
    /* Apply Abramowitz and Stegun 15.3.4 */
    F=R_pow(1-z,-a)*hyperg2F1(a,c-b,c,z/(z-1));
  }
  else if (z==1) {
    F=exp(lgammafn(c)+lgammafn(c-a-b)-lgammafn(c-a)-lgammafn(c-b));
  }
  else {
    while ((w/F)>FACCURACY) {
      j++;
      w*=((a+j-1)*(b+j-1)/(c+j-1))*(z/j);
      F+=w;
    }
  }
  return(F);
}

typedef struct C_int_struct
{
  void (*f)(double *x, int n, SEXP theta) ;    /* function */
    SEXP theta;  /*other args to f */
} C_int_struct, *IntStruct;

static void Cintfn(double *x, int n, void *ex)
{
  
  IntStruct IS = (IntStruct) ex;
  
  IS->f(x, n, IS->theta);
  
  return;
}

double phi1_int(double a, double b, double c, double x, double y, int div, double scale) {
  
  // use R's integrate code from QUADPACK to obtain marginal likelihood

    double offset = 0.0, lower = 0.0, upper = 1.0, epsabs, epsrel;
    double result, abserr, *work;
    int neval, ier, limit=200, lenw, last, *iwork;
    SEXP Rtheta;
    C_int_struct is;
    
    
    epsabs = R_pow(DBL_EPSILON, 0.25);
    epsrel = epsabs;
    lenw = 4 * limit;
    iwork = (int *) R_alloc((size_t) limit, sizeof(int));
    work = (double *) R_alloc((size_t) lenw, sizeof(double));
    
        
    PROTECT(Rtheta = allocVector(REALSXP, 7));
    REAL(Rtheta)[0] = a;
    REAL(Rtheta)[1] = b;
    REAL(Rtheta)[2] = c;
    REAL(Rtheta)[3] = x;
    REAL(Rtheta)[4] = y;
    REAL(Rtheta)[5] = (double) div;
    REAL(Rtheta)[6] = scale;
    
//    if (y <= 0 )  Rprintf("Error in Phi1: y  <= 0");
    
    is.f = phi1_density;
    is.theta = Rtheta;

    
    Rdqags(Cintfn, (void*)&is, &lower,&upper,&epsabs,&epsrel,&result,
           &abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);
   
     if (!R_FINITE(result)) {
         Rprintf("phi return: int %lf W=%lf div= %lf scale=%le \n", 
                 log(result), x, (double) div, scale);
     }
    
    if (scale > 0) offset = (double) div * log(scale);
    
    UNPROTECT(1);
    return(log(result)- offset);
}

void phi1_density(double *u, int n, SEXP Rtheta) {

  double a,b,c,x,y,z,div, scale;
  int i, j;

  PROTECT(Rtheta);
  SEXP Rex = PROTECT(duplicate(Rtheta));

  a = REAL(Rex)[0];
  b = REAL(Rex)[1];
  c = REAL(Rex)[2];
  x = REAL(Rex)[3];
  y = REAL(Rex)[4];
  div = REAL(Rex)[5]; 
  scale = REAL(Rex)[6];
  


for (i=0; i < n; i++) {
  z = u[i];
  u[i] =  exp((a - 1.0)*log(z) + (c - a - 1.0)*log(1.0 - z) - b*log(1.0 - y*z));
  for (j=0; j < (int) div; j++) {
    u[i] *= scale*exp(z*x/div);
  }
  if (!R_FINITE(u[i])) {
    Rprintf("integrate: z= %lf phi1=%lf W=%lf a=%lf b=%lf c=%lf y=%lf scale=%le div=%lf\n", z, u[i], x, a, b, c, y, scale, div);
    }
  
  u[i] *=  exp(lgammafn(c) -lgammafn(a) - lgammafn(c - a));
  }
  UNPROTECT(2);
}

double tcch_int(double a, double b, double r, double s, double v,  double k) {
  
  
  
  // use R's integrate code from QUADPACK to obtain marginal likelihood
  
  double offset = 0.0, lower = 0.0, upper,  epsabs, epsrel;
  double result, abserr, *work;
  int neval, ier, limit=200, lenw, last, *iwork;
  SEXP Rtheta;
  C_int_struct is;
  
  
  epsabs = R_pow(DBL_EPSILON, 0.25);
  epsrel = epsabs;
  lenw = 4 * limit;
  iwork = (int *) R_alloc((size_t) limit, sizeof(int));
  work = (double *) R_alloc((size_t) lenw, sizeof(double));
  
  
  PROTECT(Rtheta = allocVector(REALSXP, 8));
  REAL(Rtheta)[0] = a;
  REAL(Rtheta)[1] = b;
  REAL(Rtheta)[2] = r;
  REAL(Rtheta)[3] = s;
  REAL(Rtheta)[4] = v;
  REAL(Rtheta)[5] = k;

  
 // if (v <= 0 | v > 1)  Rprintf("Error in tcch: v is outside [0, 1)");
  upper = 1.0/v;
  
  is.f = tcch_density;
  is.theta = Rtheta;
  
  
  Rdqags(Cintfn, (void*)&is, &lower,&upper,&epsabs,&epsrel,&result,
         &abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);
  
  if (!R_FINITE(result)) {
    Rprintf("ttch return Inf: int %lf s=%lf a=%lf b=%lf r=%lf  v= %lf k=%lf lower=%lf upper=%lf\n", 
            log(result), s, a, b, r, v, k, lower, upper);
  }
  UNPROTECT(1);
  return(log(result));
}

void tcch_density(double *u, int n, SEXP Rtheta) {
  
  double a, b, r, s, v, k, z;
  int i, j;
  
  PROTECT(Rtheta);
  SEXP Rex = PROTECT(duplicate(Rtheta));
  
  a = REAL(Rex)[0];
  b = REAL(Rex)[1];
  r = REAL(Rex)[2];
  s = REAL(Rex)[3];
  v = REAL(Rex)[4];
  k = REAL(Rex)[5];
 
  
  for (i=0; i < n; i++) {
    z = u[i];
    u[i] =  (a - 1.0)*log(z) + 
            (b - 1.0)*log(1.0 - v*z) - r*log(k + (1.0 - k)*v*z);
//    for (j=0; j < (int) div; j++) {
    u[i] += -z*s;
//    }
    u[i] = exp(u[i]);
 /*   if (!R_FINITE(u[i])) {
      Rprintf("integrate: z= %lf tcch=%lf W=%lf a=%lf b=%lf r=%lf v=%lf  k = %lf, scale=%le div=%lf\n", 
              z, u[i], a, b, r, v, k, scale, div);
    }
*/
  }
  UNPROTECT(2);
}


/*--------------------------------------------------------------------------
  HyperTwo
  Function for Phi1(a,b,c,x,y)
  Assumes: 0<a<c, 0<b, y<1
  Use rule T5 if y<0.  Then use rule T1 if x>=0 or rule T2 if x<0.
--------------------------------------------------------------------------*/
double HyperTwo(double a, double b, double c, double x, double y) {
  int m=0;
  double F, zf=1.0, zfg;

  if (y<0)
    F=exp(x)*R_pow(1-y,-b)*HyperTwo(c-a,b,c,-x,y/(y-1));
  else {
    zfg=hyperg2F1(b,a,c,y);
    F=zfg;
    if (x<0) {
      while ((zfg/F)>FACCURACY) {
        m++;
        zf*=((c-a+m-1)/(c+m-1))*(-x/m);
        zfg=zf*hyperg2F1(b,a,c+m,y);
        F+=zfg;
      }
      F*=exp(x);
    }
    else {
      while ((zfg/F)>FACCURACY) {
        m++;
        zf*=((a+m-1)/(c+m-1))*(x/m);
        zfg=zf*hyperg2F1(b,a+m,c+m,y);
        F+=zfg;
      }
    }
  }
  return(F);
}


void phi1(double *a, double *b, double *c, double *x, double *y, int *div, double *scale, double *phi, 
          int *npara)
{
  int k;
  
  if (*div > 1) {
    for (k = 0; k < *npara; k++) {
      phi[k] = phi1_int(a[k], b[k], c[k], x[k], y[k], *div, *scale); 
    }
  }
  else{ 
    for (k = 0; k < *npara; k++) {
    //   if (x[k] <0) {
      /*  Since Linux system tends to report error for negative x, we
	  use the following fomular to convert it to positive value
	  1F1(a, b, x) = 1F1(b - a, b, -x) * exp(x) */
    /*      a[k] = b[k] - a[k];
      y[k] = hyperg(a[k], b[k], -x[k])*exp(x[k]);
    }
    else {
      y[k] = hyperg(a[k], b[k], x[k]);
    }*/
    phi[k] = log(HyperTwo(a[k], b[k], c[k], x[k], y[k]));
    }
  }
}

void tcch(double *a, double *b, double *r, double *s, double *v, double *theta, 
          double *tcch, 
          int *npara)
{
  int k;
  
    for (k = 0; k < *npara; k++) {
      tcch[k] = tcch_int(a[k], b[k], r[k], s[k], v[k], theta[k]); 
    }
}
