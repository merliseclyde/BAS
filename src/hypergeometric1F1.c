#include <math.h>
#include <Rmath.h>
#include <Rinternals.h>

extern double hyperg(double, double, double), lgammafn(double);
double loghyperg1F1_laplace(double, double, double);

double loghyperg1F1(double a, double b, double x, int laplace)
{
  double y;

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
    /* Laplace approximation assumes -x  for x positive
    if ( x <= 0.0 ){y = loghyperg1F1_laplace(a, b, -x); }
    else { y = loghyperg1F1_laplace(b - a, b, x) + x;} */
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
  double  mode,mode1, mode2, lprec, prec, logy;

  /* int u^(a-1) (1-u)^(b-1) exp(-x u) du   assuming that x >= 0 */

  prec = 0.0;
  logy = 0.0;

  if ( x <= 0.0) {
     if (x < 0.0) {
         x = -x;
       logy =  -lgammafn(b) - lgammafn(a) + lgammafn(a+b);

	//	mode = (2.0 - 2.0* a + b - x - sqrt(pow(b, 2.0) - 2.0*b*x + x*(4.0*(a-1.0)+x)))/
	//  (2*(a - 1.0 - b));
       mode1 = .5*(-a + b + x - sqrt( 4.*a*b + pow(a - b - x, 2.0)))/a;
       mode1 =  1.0/(1.0 + mode1);
       mode2 = .5*(-a + b + x + sqrt( 4.*a*b + pow(a - b -x, 2.0)))/a;
       mode2 =  1.0/(1.0 + mode2);
       if (a*log(mode1) + b*log(1.0 - mode1) - x*mode1  >
	   a*log(mode2) + b*log(1.0 - mode2) - x*mode2) mode = mode1;
       else mode = mode2;
       //       Rprintf("mode 1 %lf, mode %lf\n", mode1, mode2);
	if (mode < 0) {
	  mode = 0.0;
	  warning("1F1 Laplace approximation on boundary\n");
	}
        else{
	  /*	  prec = a*mode*(1.0 - mode) + (1.0-mode)*(1.0 - mode)*b +
	         x*pow(1.0-mode, 3.0) - x*mode*(1.0 -
	         mode)*(1.0-mode); */
	  prec = (1.0-mode)*((a + b - x)*pow(mode,2) + (1.0-mode)*mode*(a + b  + x));
	  if (prec > 0)  {
	      lprec = log(prec);
	      logy += a*log(mode) + b*log(1.0 - mode) - x*mode;
	      logy += -0.5*lprec + M_LN_SQRT_2PI;
	  }
	  else {prec = 0.0;}
	}

	//	Rprintf("mode %lf prec %lf, Lap 1F1(%lf, %lf, %lf) = %lf\n", mode, prec, a,b,x, logy);
     }
     else {logy = 0.0;}
  }
  else {
    logy = x + loghyperg1F1_laplace(b - a, a, -x);
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


