#include "gglm.h"

/************* Under the robust prior ***************** 
   default choice of hyper parameters:
   a = 1/2, b = 1, rho = 1 / (1 + p_M)
   
   Input:
   z: MCMC samples of z = 1 / (1 + g), length = iter
   invIbeta: inverse of information of beta based on all data
   betamle: mle of beta
 
   Output:
   beta: matrix, dim = c(p, iter)
*/

void robust_postbeta(int *p, double *invIbeta, double *betamle, 
    int *iter, double *z, double *beta)
{
  int t, one = 1, inc = 1, p2 = (*p) * (*p);
  double mu[*p], sigma[(*p) * (*p)];
  
  for (t = 0; t < *iter ; t++) {
    // mu = z[t] * betamle
    memcpy(mu, betamle, (*p) * sizeof(double));
    F77_CALL(dscal)(p, &(z[t]), mu, &inc);
    // sigma = z[t] * invIbeta
    memcpy(sigma, invIbeta, (*p) * (*p) * sizeof(double));
    F77_CALL(dscal)(&p2, &(z[t]), sigma, &inc);
    // beta ~ N(mu, sigma)
    r_multi_norm(p, &one, mu, sigma, &(beta[(*p) * t]));
  }
  
}

 
