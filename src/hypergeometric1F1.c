#include <math.h>
#include <Rmath.h>

extern double hyperg(double, double, double);
double hyperg1F1_laplace(double, double, double);
double hyperg1F1_laplace_pos(double, double, double);

double hyperg1F1(double a, double b, double x)
{ 
  double y, ly;


    if (x <0) {
      /*  Since Linux system tends to report error for negative x, we
	  use the following fomular to convert it to positive value 
	  1F1(a, b, x) = 1F1(b - a, b, -x) * exp(x) */
      y = hyperg(b-a, b, -x)*exp(x);
    }
    else {
      y = hyperg(a, b, x);
    }

    if (x >= 0) {
    ly = hyperg1F1_laplace(a,b,-x);
}
    else {
      ly = hyperg1F1_laplace(a,b,-x);
    }
    Rprintf("LOG 1F1(%lf, %lf, %lf) = %lf (%lf)\n", a,b,x,log(y), y);

    ly = hyperg1F1_laplace_pos(a,b,x);
    Rprintf("LOG Pos 1F1(%lf, %lf, %lf) = %lf (%lf)\n", a,b,x,log(y), y);

    return(y);
} 

double hyperg1F1_laplace(double a, double b, double x)
{ 
  double y, mode, mode_neg, lprec, logy;

  
  mode = (-2.0 + a - b - x + sqrt(pow(a, 2.0) + pow(b, 2.0) + 2.0*a*(b-x) + 2.0*b*x + x*(4+x)))/
    (2*(1.0 + b));
  
  
  Rprintf("Laplace Approx: Mode %lf \n", mode);
  if (mode < 0) {
    Rprintf("uh oh mode is neg\n");
  }

  lprec =  -((1.0 - a)/pow(mode,2.0) - 2.0*x*mode/pow((1.0 + mode),3.0) -
	     (a + b - 2.0*x)/pow((1.0 + mode),2.0));
  Rprintf("prec %lf\n", lprec);
  if (lprec > 0) lprec = log(lprec);
  else Rprintf("wrong root\n");
  
	  
  logy = (a - 1.0)*log(mode) - (a + b)*log(1.0 + mode) - x*mode/(1.0 + mode);
  logy += -0.5*lprec + M_LN_SQRT_2PI - lbeta(a, b);
//  y = hyperg1F1(a, b, x);
  
  Rprintf("LOG Laplace 1F1(%lf, %lf, %lf) = %lf (%lf)\n", a,b,x,logy, exp(logy));
	     return(logy);	
}

double hyperg1F1_laplace_pos(double a, double b, double x)
{ 
  double y, mode, mode_neg, lprec, logy;

  if ( x >= 0.0)
    {
      mode = (2.0 - 2.0* a + b - x + sqrt(pow(b, 2.0) - 2.0*b*x + x*(4.0*(a-1)+x)))/
	(2*(a - 1.0 - b));
  
  
      Rprintf("1F1 Lap Pos Mode1 %lf \n", mode);
      if (mode < 0) {
	Rprintf("uh oh mode is neeg\n");
      }

      lprec =  -((1.0 - a)/pow(mode,2.0) + 2.0*x*mode/pow((1.0 + mode),3.0) -
		 (b + 2.0*x)/pow((1.0 + mode),2.0));
      Rprintf("prec %lf\n", lprec);
      if (lprec > 0) lprec = log(lprec);
      else Rprintf("wrong root\n");

      logy =  lgamma(b) - lgamma(a) - lgamma(b-a);
  if (x > 0.0) {	  
    logy += (a - 1.0)*log(mode) - (a + b)*log(1.0 + mode) - x*mode/(1.0 + mode);
    logy += -0.5*lprec + M_LN_SQRT_2PI;
  }}
  else {
    logy = x + hyperg1F1_laplace_pos(b - a, b, -x);
  }
//  y = hyperg1F1(a, b, x);
  
  Rprintf("LOG Laplace Pos 1F1(%lf, %lf, %lf) = %lf (%lf)\n", a,b,x,logy, exp(logy));
	     return(logy);	
}


void hypergeometric1F1(double *a, double *b, double *x, double *y, int *npara)
{ 
  int k;
  for (k = 0; k < *npara; k++) {
    if (x[k] <0) {
      /*  Since Linex system tends to report error for negative x, we
	  use the following fomular to convert it to positive value 
	  1F1(a, b, x) = 1F1(b - a, b, -x) * exp(x) */
      a[k] = b[k] - a[k];
      y[k] = hyperg(a[k], b[k], -x[k])*exp(x[k]);
    }
    else {
      y[k] = hyperg(a[k], b[k], x[k]);
    }
  }
}

double shrinkage_chg(double a, double b, double s, double Q) {

  double shrinkage;
  /* Beta(a/2,(b+2)/2) 1F1(a/2,(b+2)/2,(s+Q)/2 /
     Beta(a/2,b/2) 1F1(a/2,b/2,(s+Q)/2
		       */		       
  shrinkage = exp( lbeta(a/2.0, (b+2.0)/2.0) +
		   log(hyperg(a/2.0, b/2.0 + 1.0, (s + Q)/2.0)) -
		   lbeta(a/2.0, (b+2.0)/2.0) -
		   log(hyperg(a/2.0, b/2.0, (s + Q)/2.0)));	

  //  Rprintf("shrinkage_chg:  %lf\n", shrinkage);
  if (shrinkage > 1.0)  shrinkage = 1.0;
  else if (shrinkage < 0.0) shrinkage = 0.0;
  
  return (shrinkage);
}
