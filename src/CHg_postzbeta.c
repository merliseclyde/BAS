#include "gglm.h"

/************* Under CH-g prior ******************** 
   z is a matrix of posterior samples of z = g / (1 + g),
   which has distribution CH( a / 2, (b + p) / 2, -(s + Q) / 2 );
   and its length is iter.
   beta is an array, each column is a posterior sample of coefficient
   
   Input:
   a, b, s can be vectors with the same length = npara
 
   Output:
   z: matrix, nrow = iter, ncol = npara
   beta: array, dim = c(p, iter, npara)
*/

void CHg_postzbeta(double *a, double *b, double *s, int *npara, int *p, double *Q, 
    double *invIbeta, double *betamle, int *iter, double *z, double *beta)
{
  int t, k, n = 1, inc = 1, p2 = (*p) * (*p);
  double z1, z2, temp, loga, logu, dp = (double)(*p), mu[*p], sigma[(*p) * (*p)];
  
  for(k=0; k < *npara; k++){
    
    // initial value of z
    z[k * (*iter)] = rbeta(a[k] / 2.0, (b[k] + dp) / 2.0);  
    // initial value mu = z[0] * betamle
    memcpy(mu, betamle, (*p) * sizeof(double));
    F77_CALL(dscal)(p, &(z[k * (*iter)]), mu, &inc);
    // initial value sigma = z[0] * invIbeta
    memcpy(sigma, invIbeta, (*p) * (*p) * sizeof(double));
    F77_CALL(dscal)(&p2, &(z[k * (*iter)]), sigma, &inc);
    // beta ~ N(mu, sigma)
    r_multi_norm(p, &n, mu, sigma, &(beta[k * (*iter) * (*p)]));
  
    // Metropolis-Hastings: from (t) to (t+1)
    for (t = 0; t < (*iter -1) ; t++) {
      z1 = z[t + k * (*iter)]; // curent z
      temp = rnorm( log(z1 / (1.0 - z1)), 0.5 );
      z2 = 1 / (1 + exp(-temp));
      
      loga = (a[k] / 2.0) * (log(z2) - log(z1)) + (b[k] + dp) / 2.0 * (log(1.0 - z2) - log(1.0 - z1)) + (s[k] + *Q) / 2 * (z2 - z1);
      logu = log(runif(0.0, 1.0));
    
      if(loga > logu){
        z[t + 1 + k * (*iter)] = z2;
      }
      if(loga <= logu){
        z[t + 1 + k * (*iter)] = z1;
      }
    
      // mu = z[t+1] * betamle
      memcpy(mu, betamle, (*p) * sizeof(double));
      F77_CALL(dscal)(p, &(z[t + 1 + k * (*iter)]), mu, &inc);
      // sigma = z[t+1] * invIbeta
      memcpy(sigma, invIbeta, (*p) * (*p) * sizeof(double));
      F77_CALL(dscal)(&p2, &(z[t + 1 + k * (*iter)]), sigma, &inc);
      // beta ~ N(mu, sigma)
      r_multi_norm(p, &n, mu, sigma, &(beta[(*p) * (t + 1) + k * (*iter) * (*p)]));
    }
    
  }
  
}

 
