//#include "gglm.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Linpack.h>

/* p: dimension of each sample
   n: number of samples to generate
   Mu: mean vector, length = p
   Sigma: covaraince vector, p by p non-singular symmetric matrix
   x: p by n, each column is a sample
*/
/* for matrix calculation */
static const char  *TRANSNO = "N", *UPLO = "U";
static const double ONE = 1.0, ZERO = 0.0;

void r_multi_norm(int *p, int *n, double *Mu, double *Sigma, double *x)
{
  int i, j, lwork, info = 1, inc = 1;
  double  work, tmp_ev,
      tmp_sigma[(*p) * (*p)], mu[*p],
      ev[*p], sigma[(*p) * (*p)], indptnorm[*p];
  char *jobz = "V";

  for (i=0; i< *n; i++) {

    memcpy(sigma, Sigma, (*p) * (*p) * sizeof(double));
    memcpy(mu, Mu, (*p) * sizeof(double));

    // calculate sqrt(sigma), using eigenvalue deposition:
    // sigma = U * Gamma * U^T, and columns of U are eigen vectors
    lwork = -1;
    F77_CALL(dsyev)(jobz, UPLO, p, sigma, p, ev, &work, &lwork, &info);
    lwork = (int)(work);
    double work2[lwork];
    F77_CALL(dsyev)(jobz, UPLO, p, sigma, p, ev, work2, &lwork, &info);


    for (j=0; j< *p; j++) {
      tmp_ev = pow(ev[j], 0.25);
      F77_CALL(dscal)(p, &tmp_ev, &sigma[j * (*p)], &inc);
      //generate indptnorm, a vector of p independent standard normals
      indptnorm[j] = rnorm(0.0, 1.0);
    }

    //each column of x = sqrt(Sigma) %*% indptnorm + mu = sigma %*% sqrt(ev) %*% t(sigma) %*% indeptnorm + mu
    F77_CALL(dsyrk)(UPLO, TRANSNO, p, p, &ONE, sigma, p, &ZERO, tmp_sigma, p);
    F77_CALL(dsymv)(UPLO, p, &ONE, tmp_sigma, p, indptnorm, &inc, &ONE, mu, &inc);

    memcpy(&x[i * (*p)], mu, (*p)*sizeof(double));
  }

}


