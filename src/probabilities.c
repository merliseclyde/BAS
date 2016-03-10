#include "sampling.h"

void compute_modelprobs(SEXP Rmodelprobs,  SEXP Rlogmarg, SEXP Rpriorprobs, int k)
{
	int m;
	double nc, bestmarg, *modelprobs, *logmarg, *priorprobs;

	logmarg = REAL(Rlogmarg);
	modelprobs = REAL(Rmodelprobs);
	priorprobs = REAL(Rpriorprobs);

	bestmarg = logmarg[0];
	nc = 0.0;

	for (m = 0; m < k; m++) {
		if (logmarg[m] > bestmarg) bestmarg = logmarg[m];
	}

	for (m = 0; m < k; m++) {
		modelprobs[m] = logmarg[m] - bestmarg;
		nc += exp(modelprobs[m])*priorprobs[m]; 
	}	

	for (m = 0; m < k; m++) {
	  /*		modelprobs[m] = exp(modelprobs[m] +
			log(priorprobs[m]) - log(nc)); */
		modelprobs[m] = exp(modelprobs[m] - log(nc))*priorprobs[m];
	}	
}


void compute_margprobs(SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, double *margprobs, int k, int p)
{
	int m, j, *model;
	double *modelprobs;
	modelprobs = REAL(Rmodelprobs);
	for (j=0; j< p; j++)  margprobs[j] = 0.0;
	for(m=0; m< k; m++) {
		model = INTEGER(VECTOR_ELT(modelspace,m));
		for (j = 0; j < INTEGER(modeldim)[m]; j ++) {
			margprobs[model[j]] += modelprobs[m];
		} 
	}
}



void compute_margprobs_old(Bit **models, SEXP Rmodelprobs, double *margprobs, int k, int p)
{
  int m, j;
  double *modelprobs;
  modelprobs = REAL(Rmodelprobs);

  for (j=0; j< p; j++) {
    margprobs[j] = 0.0;
   for(m=0; m< k; m++) {
     if (models[m][j])
        margprobs[j] += modelprobs[m];
      }
    }
}

double compute_prior_probs(int *model, int modeldim, int p, SEXP modelprior) {
  const char *family;
  double *hyper_parameters, priorprob = 1.0;


  family = CHAR(STRING_ELT(getListElement(modelprior, "family"),0));
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));


  if  (strcmp(family, "Beta-Binomial") == 0) 
    priorprob = beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Beta-Binomial") == 0) 
    priorprob = trunc_beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Poisson") == 0) 
    priorprob = trunc_poisson(modeldim, p, hyper_parameters);
  if (strcmp(family, "Bernoulli") == 0) 
    priorprob = Bernoulli(model, p, hyper_parameters);
  return(priorprob);
}
     
double Bernoulli(int *model, int p, double *hyper) {
  double prior;
  int j;

  for (j=1, prior=1.; j < p; j++) {
    switch(model[j]) {
    case 0:
      prior *= (1. - hyper[j]);
      break;
    case 1:
      prior *= hyper[j];
      break;
    default:
      prior *= 1.;
      break;
      }
  }
  return(prior);
}


double beta_binomial(int modeldim, int p, double *hyper) {
  /* modeldim and p include the intercept so subtact 1 from each */ 
  return(exp(lbeta((double) modeldim - 1.0 + hyper[0], (double) (p - modeldim) + hyper[1]) -
	     lbeta(hyper[0], hyper[1])));
}

double trunc_beta_binomial(int modeldim, int p, double *hyper) {
  /* modeldim and p include the intercept so subtact 1 from each */

  double prior;
  if ((double) (modeldim -1) <= hyper[2]) {
      prior = exp(lbeta((double) modeldim - 1.0 + hyper[0], (double) (p - modeldim) + hyper[1]) - 
		  lbeta(hyper[0], hyper[1]));
      //      Rprintf("pass \n");
      
    }
  else {prior = 0.0;}

  /*  Rprintf("prior %lf pmodel= %d s0 = %lf\n", prior, modeldim,
      hyper[2]); */
  return(prior);
}

double trunc_poisson(int modeldim, int p, double *hyper) {
  /* modeldim and p include the intercept so subtact 1 from each */

  double prior = 0.0;
  if ((double) (modeldim -1) <= hyper[1]) {
      prior = dpois(modeldim - 1, hyper[0], 0);
    }
  return(prior);
}

