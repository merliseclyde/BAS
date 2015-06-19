#include <math.h>
extern double hyperg(double, double, double);

double hyperg1F1(double a, double b, double x)
{ 
  double y;
    if (x <0) {
      /*  Since Linex system tends to report error for negative x, we
	  use the following fomular to convert it to positive value 
	  1F1(a, b, x) = 1F1(b - a, b, -x) * exp(x) */
      y = hyperg(b-a, b, -x)*exp(x);
    }
    else {
      y = hyperg(a, b, x);
    }
    //    Rprintf("LOG 1F1(%lf, %lf, %lf) = %lf (%lf)\n", a,b,x,log(y), y);
    return(y);
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

