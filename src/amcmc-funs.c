// AMC Helper functions

#include "bas.h"

void  update_Cov(double *Cov, double *priorCov, double *SSgam, double *marg_probs, int n, int m, int print) {
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
  // Compute SSgam - m prob prob^T 
  F77_NAME(dsyr)("U", &n,  &alpha, &marg_probs[0], &inc,  &Cov[0], &n FCONE);
  alpha = 1.0/(double) m;
  // divide SSgam - m prob prob^T by m -> Cov
  F77_NAME(dscal)(&matsize,  &alpha,  &Cov[0], &inc);
  // Add priorCov to Cov
  F77_NAME(daxpy)(&matsize, &one, &priorCov[0], &inc, &Cov[0], &inc);
  
  if (print == 1) {
    Rprintf("Cov:\n");
    
    for (j=0; j < n; j++ ){
      for(i = 0; i < n; i++){
        Rprintf("%f ", Cov[i*n + j]);}
      Rprintf("\n");
    }
  }
  F77_NAME(dpotrf)("U", &n,  &Cov[0], &n, &info FCONE);
  F77_NAME(dtrtri)("U","N", &n, &Cov[0], &n, &info FCONE FCONE);
  
  if (print == 1) {
    Rprintf("inverse of Chol(Cov(SSgam)):\n");
    for (j=0; j < n; j++ ){
      for(i = 0; i < n; i++){
        Rprintf("%f ", Cov[i*n + j]);}
      Rprintf("\n");
    } 
  }
}


double cond_prob(double *model, int j, int n, double *mean, double *beta_matrix , double delta) {
  double  prob;
  int i;
  
  prob = mean[j]; 
  //  Rprintf("%f\n", prob);
  for (i = 0; i < j; i++) {
  // old code had a - beta ??? should be +
    prob += beta_matrix[j*n + i]*(model[i] - mean[i])/beta_matrix[j*n+j];
    //    Rprintf("model %f beta %f mean %f ", model[i], beta_matrix[i*n + j], mean[i]);
  }
  //  Rprintf("\n%f ", prob);
  if (prob <= 0.0) prob = delta;
  if (prob >= 1.0) prob = 1.0 -  delta;
  //    Rprintf("%f \n", prob);
  return(prob);
}
