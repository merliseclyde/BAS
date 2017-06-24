#include <stdio.h>
#include <stdlib.h>
#include "mconf.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


/* for matirx calculation */
static const char *TRANSYES = "T", *TRANSNO = "N", *UPLO = "U";
static const double ONE = 1.0, ZERO = 0.0;


////////// functions
double hyperg( double a, double b, double x);
void r_multi_norm(int *p, int *n, double *Mu, double *Sigma, double *x);

//void postzbeta(double *a, double *b, double *s, int *npara, int *p, double *Q,
//               double *invIbeta, double *betamle, int *iter, double *z, double *beta);
