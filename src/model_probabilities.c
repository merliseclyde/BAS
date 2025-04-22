// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include "bas.h"

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

void compute_modelprobs_HT(SEXP Rmodelprobs,  SEXP Rlogmarg, SEXP Rpriorprobs, 
                           SEXP Rsampleprobs, int k, int MCsamples)
{
  int m;
  double nc, bestmarg, *modelprobs, *logmarg, *priorprobs, *sampleprobs;
  
  logmarg = REAL(Rlogmarg);
  modelprobs = REAL(Rmodelprobs);
  priorprobs = REAL(Rpriorprobs);
  sampleprobs = REAL(Rsampleprobs); 
  bestmarg = logmarg[0];
  nc = 0.0;
  
  for (m = 0; m < k; m++) {
    if (logmarg[m] > bestmarg) bestmarg = logmarg[m];
    if (sampleprobs[m] > 0.0) modelprobs[m]  =  -log(1.0 - (pow(1.0 - sampleprobs[m], (double) MCsamples)));
  }
  
  for (m = 0; m < k; m++) {
    if (sampleprobs[m] > 0.0) {
      modelprobs[m] += logmarg[m] - bestmarg;
      nc += exp(modelprobs[m])*priorprobs[m];
    }
  }
  
  for (m = 0; m < k; m++) {
    if (sampleprobs[m] > 0.0) modelprobs[m] = exp(modelprobs[m] - log(nc))*priorprobs[m];
    else {modelprobs[m] = 0.0;}
  }
}

void compute_modelprobs_Bayes_HT(SEXP Rmodelprobs,  SEXP Rlogmarg, SEXP Rpriorprobs, 
                           SEXP Rsampleprobs, int M)
{
  int m;
  double nc, bestmarg, *modelprobs, *logmarg, *priorprobs, *sampleprobs, 
         HT = 0.0, probinS = 0.0;
  
  logmarg = REAL(Rlogmarg);
  modelprobs = REAL(Rmodelprobs);
  priorprobs = REAL(Rpriorprobs);
  sampleprobs = REAL(Rsampleprobs); 
  bestmarg = logmarg[0];
  nc = 0.0;
  
  for (m = 0; m < M; m++) {
    if (logmarg[m] > bestmarg) bestmarg = logmarg[m];
  }
  
  for (m = 0; m < M; m++) {
      modelprobs[m] = logmarg[m] - bestmarg + log(priorprobs[m]);
      if (sampleprobs[m] > 0.0) {
        HT += exp(modelprobs[m] - log(sampleprobs[m]));  }
      probinS += sampleprobs[m];
      nc += exp(modelprobs[m]);
  }
  
  nc += (1.0 - probinS)*HT/ (double) M;
  for (m = 0; m < M; m++) {
    modelprobs[m] = exp(modelprobs[m] - log(nc));
  }
}

void compute_margprobs(SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, double *margprobs, 
                       int k, int p)
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

void compute_margprobs_Bayes_BAS_MCMC(SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, double *margprobs, SEXP Rsampleprobs, 
                       int M, int p)
{
  int m, j, *model;
  double *modelprobs;
  double *correction_probs, *probNotInS_j;
  
  correction_probs = (double *) R_alloc(p, sizeof(double));;
  probNotInS_j = (double *) R_alloc(p, sizeof(double));
  memset(correction_probs, 0, p * sizeof(double));
  memset(probNotInS_j, 0, p * sizeof(double));
  modelprobs = REAL(Rmodelprobs);
  
  for (j=0; j< p; j++)  {
    probNotInS_j[j] = margprobs[j];  // on entry margprobs[j] is the marginal sampling probability
    margprobs[j] = 0.0;
    
  }
  for(m=0; m< M; m++) {
    model = INTEGER(VECTOR_ELT(modelspace,m));
    
    for (j = 0; j < INTEGER(modeldim)[m]; j ++) {
      margprobs[model[j]] += modelprobs[m];
      probNotInS_j[model[j]] -= REAL(Rsampleprobs)[m];
      correction_probs[model[j]] += modelprobs[m]/REAL(Rsampleprobs)[m];
    }
  }
  for (j=0; j< p; j++)  margprobs[j] += correction_probs[j]*(probNotInS_j[j])/((double) M);
}

void compute_sampleprobs_modelspace_Bernoulli(SEXP modelspace, SEXP modeldim, SEXP Rsampleprobs, SEXP Rprobs, 
                       int nModels, int p)
{
  int m, j, *model; 
  int *modelVec;
  modelVec = ivecalloc(p);
  memset(modelVec, 0, p * sizeof(int));
  
    for(m=0; m < nModels; m++) {
    memset(modelVec, 0, p * sizeof(int));
    model = INTEGER(VECTOR_ELT(modelspace,m));
    for (j = 0; j < INTEGER(modeldim)[m]; j ++) {
      modelVec[model[j]] = 1.0;
    }
    REAL(Rsampleprobs)[m] = compute_sample_probs_bernoulli(Rprobs, modelVec, p);

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

int no_prior_inclusion_is_1(int p, double *probs) {

  int noInclusionIs1 = 0;
  // loop starts from 1 since the intercept is corrected for in the model prior functions
  for (int i = 1; i < p; i++) { 
  	if (probs[i] > (1.0 - DBL_EPSILON)) {
  		noInclusionIs1++;
  	}
  }
  return noInclusionIs1;
}

double compute_sample_probs_bernoulli(SEXP Rprobs, int *model, int p) {
  int j;
  double pigamma = 1.0;
  for (j = 0; j < p; j++) {
    pigamma *= ((double) model[j])*REAL(Rprobs)[j] + (1.0 - ((double) model[j]))*(1.0 -  REAL(Rprobs)[j]);
  }

return(pigamma);
}

double compute_prior_probs(int *model, int modeldim, int p, SEXP modelprior, int noInclusionIs1) {
  const char *family;
  double *hyper_parameters, priorprob = 1.0;


  family = CHAR(STRING_ELT(getListElement(modelprior, "family"),0));
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));

  // do not reduce p by the number of predictors that are always included
  // Gitub issue # 87
  if (strcmp(family, "Bernoulli") == 0)
    priorprob = Bernoulli(model, p, hyper_parameters);
  
  // reduce the model space by the number of predictors that are always included 
  p -= noInclusionIs1;
  modeldim -= noInclusionIs1;

  if  (strcmp(family, "Beta-Binomial") == 0)
    priorprob = beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Beta-Binomial") == 0)
    priorprob = trunc_beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Poisson") == 0)
    priorprob = trunc_poisson(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Power-Prior") == 0)
    priorprob = trunc_power_prior(modeldim, p, hyper_parameters);
// Need to add
//  if (strcmp(family, "Hereditary") == 0)
//    priorprob = Hereditary(model, p, hyper_parameters);
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
/*  Can't ever get here
    default:
      prior *= 1.;
      break; 
 */
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
  /* modeldim and p include the intercept so subtract 1 from each */

  double prior = 0.0;
  if ((double) (modeldim -1) <= hyper[1]) {
      prior = dpois(modeldim - 1, hyper[0], 0);
    }
  else {prior = 0.0;}

  return(prior);
}

double trunc_power_prior(int modeldim, int p, double *hyper) {
  /* modeldim and p include the intercept so subtract 1 from each */

  double prior = 0.0;
  if ((double) (modeldim -1) <= hyper[1]) {
    prior = exp(-((double) modeldim - 1.0)*((double) hyper[0])*log((double) p));
    }
  else {prior = 0.0;}

  return(prior);
}

