/**********************************************************************************
 * This file is part of the R package: G-prior for Generalized Linear Models (gglm)
 **********************************************************************************  
 * gglm is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 2 of the License, or any later version.
 *
 * gglm is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * Main function: to be filled
 * Goal: to be filled
 *
 * version  12/24/2012 
 * Questions? Contact Yingbo Li (yl118@duke.edu)
 *
 ********************************************************************************/

/*
#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
*/

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
void hypergeometric1F1(double *a, double *b, double *x, double *y, int *npara);
void r_multi_norm(int *p, int *n, double *Mu, double *Sigma, double *x);

//void postzbeta(double *a, double *b, double *s, int *npara, int *p, double *Q, 
//               double *invIbeta, double *betamle, int *iter, double *z, double *beta);