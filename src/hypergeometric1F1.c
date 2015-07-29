#include <math.h>
#include <Rmath.h>
#include <Rinternals.h>

extern double hyperg(double, double, double);
double loghyperg1F1_laplace(double, double, double);


double loghyperg1F1(double a, double b, double x, int laplace)
{ 
  double y, ly;

  if (laplace == 0) {
    if (x <0) {
      /*  Since Linux system tends to report error for negative x, we
	  use the following fomular to convert it to positive value 
	  1F1(a, b, x) = 1F1(b - a, b, -x) * exp(x) */
      y = log(hyperg(b-a, b, -x)) + x;
    }
    else {
      y = log( hyperg(a, b, x));
    }
  }
  else {
    y = loghyperg1F1_laplace(a, b, x);
  }
    
    //    Rprintf("LOG Cephes 1F1(%lf, %lf, %lf) = %lf (%lf)\n", a,b,x,log(y), y);
    

    //    ly = hyperg1F1_laplace(a,b,x);
    //    Rprintf("called from hyperg1F1: LOG Pos 1F1(%lf, %lf, %lf) = %lf (%lf)\n", a,b,x,ly, exp(ly));
  if (!R_FINITE(y)  &&  laplace == 0) {
    warning("Cephes 1F1 function returned NA, using Laplace approximation");
    y = loghyperg1F1_laplace(a, b, x);  // try Laplace approximation
  }
  
  return(y);
} 


double loghyperg1F1_laplace(double a, double b, double x)
{ 
  double  mode, lprec, prec, logy;

  logy =  lgamma(b) - lgamma(a) - lgamma(b-a);

  if ( x >= 0.0)
    {
      
      if (x > 0.0) {
 	mode = (2.0 - 2.0* a + b - x - sqrt(pow(b, 2.0) - 2.0*b*x + x*(4.0*(a-1.0)+x)))/
	  (2*(a - 1.0 - b));
  
  
	if (mode < 0) {
	  mode = 0.0;
	  Rprintf("1F1 Laplace approximation on boundary\n");
	}

	prec =  -((1.0 - a)/pow(mode,2.0) + 2.0*x*mode/pow((1.0 + mode),3.0) +
		 (b - 2.0*x)/pow((1.0 + mode),2.0));
	if (prec > 0) lprec = log(prec);
	else prec = 0.0;
      
      
	logy += (a - 1.0)*log(mode) - b*log(1.0 + mode) + x*mode/(1.0 + mode);
	logy += -0.5*lprec + M_LN_SQRT_2PI;
	logy += pnorm(0.0, mode, 1.0/sqrt(prec), 0, 1);

	//	Rprintf("mode %lf prec %lf, Lap 1F1(%lf, %lf, %lf) = %lf", mode, prec, a,b,x); 
      }}
  else {
    logy = x + loghyperg1F1_laplace(b - a, b, -x);
  }

  return(logy);	
}


void hypergeometric1F1(double *a, double *b, double *x, double *y, int *npara, int *Method)
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
    y[k] = loghyperg1F1(a[k], b[k], x[k], Method[k]);
  }   
}

double shrinkage_chg(double a, double b, double Q, int laplace) {

  double shrinkage;
  /* Beta(a/2,(b+2)/2) 1F1(a/2,(b+2)/2,(s+Q)/2 /
     Beta(a/2,b/2) 1F1(a/2,b/2,(s+Q)/2
  */		       
  /* shrinkage = exp( lbeta(a/2.0, (b+2.0)/2.0) +
		   log(hyperg1F1(a/2.0, b/2.0 + 1.0, Q/2.0)) -
		   lbeta(a/2.0, (b)/2.0) -
		   log(hyperg1F1(a/2.0, b/2.0, Q/2.0)));
   */
   //    Rprintf("shrinkage_chg:  %lf\n", shrinkage);
    shrinkage = exp( lbeta(a/2.0, b/2.0 + 1.0) +
		     loghyperg1F1(a/2.0, b/2.0 + 1.0,  Q/2.0, laplace) -
		     lbeta(a/2.0, b/2.0) -
		     loghyperg1F1(a/2.0, b/2.0,  Q/2.0, laplace));	

    //Rprintf("Laplace shrinkage_chg:  %lf\n", shrinkage);
  if (shrinkage > 1.0)  shrinkage = 1.0;
  else if (shrinkage < 0.0) shrinkage = 0.0;
  
  return (shrinkage);
}
