// AMC Helper functions

#include "bas.h"

void  update_Cov(double *Cov, double *priorCov, double *SSgam, double *marg_probs, double lambda, int n, int m, int print) {
  double alpha, one=1.0;
  int matsize=n*n, inc=1, i,j,info=1;
  
  memcpy(Cov, SSgam, sizeof(double)*matsize);
  if (print == 1) {
    Rprintf("SS: %d iterations\n", m);
    
    for (j=0; j < n; j++ ){
      Rprintf("%d  %f  ", j, marg_probs[j]);
      for(i = 0; i < n; i++){
        Rprintf("%f ", Cov[i*n + j]);}
      Rprintf("\n");
    }
  }
  alpha  = -(double) m;
  // Compute SSgam - m prob prob^T  = SS_n
  F77_NAME(dsyr)("U", &n,  &alpha, &marg_probs[0], &inc,  &Cov[0], &n FCONE);
  // Add priorSS to observed SS = SS_0 + SS_n
  F77_NAME(daxpy)(&matsize, &one, &priorCov[0], &inc, &Cov[0], &inc);
  
  alpha = 1.0/((double) m + lambda - (double) n - 1.0);
  // divide SS_0 + SS by (m + lambda - n - 1) to get posterior expectation of Sigma
  F77_NAME(dscal)(&matsize,  &alpha,  &Cov[0], &inc);
  
  if (print == 1) {
    Rprintf("Cov:\n");
    
    for (j=0; j < n; j++ ){
      for(i = 0; i < n; i++){
        Rprintf("%f ", Cov[i*n + j]);}
      Rprintf("\n");
    }
  }
  // compute the cholesky of Cov = U^T U
  F77_NAME(dpotrf)("U", &n,  &Cov[0], &n, &info FCONE);
  // compute the inverse of Cholesky factor U.inv. 
  F77_NAME(dtrtri)("U","N", &n, &Cov[0], &n, &info FCONE FCONE);
  // verified correct 3/15/2024
  if (print == 1) {
    Rprintf("inverse of Chol(Cov(SSgam)):\n");
    for (j=0; j < n; j++ ){
      for(i = 0; i < n; i++){
        Rprintf("%f ", Cov[i*n + j]);}
      Rprintf("\n");
    } 
  }
}


double cond_prob(double *model, int j, int n, double *mean, double *u_inv , double delta) {
  double  prob;
  int i;
  
  // betas =  I - diag(U^{-1}) U^-T  
  // Note that elements of u_inv is U^{-T} due to FORTRAN column order
  // E[Y] = betas (Ybar)
  // E[Y_j | Y_<j ]   =  Ybar_j + betas(Y<j - Ybar<j)
  prob = mean[j]; 
  //  Rprintf("%f\n", prob);
  for (i = 0; i < j; i++) {
    prob += - u_inv[j*n + i]*(model[i] - mean[i])/u_inv[j*n+j];
   // Rprintf("j %d model %f beta %f mean %f \n", j, model[i], -u_inv[j*n + i]/u_inv[j*n + j], mean[i]);
  }
  // Rprintf("\n%f ", prob);
  if (prob <= 0.0) prob = delta;
  if (prob >= 1.0) prob = 1.0 -  delta;
  //    Rprintf("%f \n", prob);
  return(prob);
}
