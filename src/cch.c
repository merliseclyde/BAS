/*--------------------------------------------------------------------------
  CONHYPTWO1(a,b,c,x,y)
  Phi1 hypergeometric function of two variables.
  Gradshteyn and Ryzhik (1965), equation 9.261.1
  This program assumes:  0<a<c, 0<b, y<=1.
  x is a column vector of real numbers, other parameters scalar.
--------------------------------------------------------------------------*/
#include <math.h>
#include <Rmath.h>
// #include "mex.h"
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


void phi1(double *a, double *b, double *c, double *x, double *y, double *phi, int *npara)
{
  int k;
  for (k = 0; k < *npara; k++) {
    //   if (x[k] <0) {
      /*  Since Linex system tends to report error for negative x, we
	  use the following fomular to convert it to positive value
	  1F1(a, b, x) = 1F1(b - a, b, -x) * exp(x) */
    /*      a[k] = b[k] - a[k];
      y[k] = hyperg(a[k], b[k], -x[k])*exp(x[k]);
    }
    else {
      y[k] = hyperg(a[k], b[k], x[k]);
    }*/
    phi[k] = HyperTwo(a[k], b[k], c[k], x[k], y[k]);
  }
}

