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
                           SEXP Rsampleprobs, int M, double *eta, double *nc)
{
  int m;
  double bestmarg, *modelprobs, *logmarg, *priorprobs, *sampleprobs, 
         HT = 0.0, probinS = 0.0, correction = 0.0, nmodels = 0.0;
  
  logmarg = REAL(Rlogmarg);
  modelprobs = REAL(Rmodelprobs);
  priorprobs = REAL(Rpriorprobs);
  sampleprobs = REAL(Rsampleprobs); 
  bestmarg = logmarg[0];

  
  for (m = 0; m < M; m++) {
    if (logmarg[m] > bestmarg) bestmarg = logmarg[m];
  }
  
  for (m = 0; m < M; m++) {
    if (sampleprobs[m] > 0.0) {
      modelprobs[m] = logmarg[m] - bestmarg + log(priorprobs[m]);
      probinS += sampleprobs[m];
      HT += exp(modelprobs[m] - log(sampleprobs[m]));  
      *nc += exp(modelprobs[m]);
      nmodels += 1.0;
    }
  }  
    
  *eta = HT/nmodels;
  correction = *eta*(1.0 - probinS);
  Rprintf("eta = %lf probinS = %lf  NC = %lf  correction = %lf", *eta, probinS, *nc, correction);
  *nc += (1.0 - probinS)* *eta;
  Rprintf(" corrected NC = %lf \n", *nc);
  
  for (m = 0; m < M; m++) {
    if (sampleprobs[m] > 0.0) {
      modelprobs[m] = exp(modelprobs[m] - log(*nc));
    }
    else {modelprobs[m] = 0.0;}
    
  }
}

void compute_margprobs(SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, double *margprobs, 
                       int k, int p) {
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

void compute_margprobs_Bayes_BAS_MCMC(SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, SEXP Rprobs, SEXP Rsampleprobs, 
                       int M, int p, double eta, double NC)
{
  int m, j, *model, *n;
  double *modelprobs;
  double *beta, probNotInS = 1.0;
  SEXP samplemargs = PROTECT(duplicate(Rprobs)); 
  modelprobs = REAL(Rmodelprobs);
  
  for (j=0; j< p; j++) {
   // Rprintf("j = %d uncorrected pip = %lf \n", j, REAL(Rprobs)[j]);
    REAL(Rprobs)[j] = 0.0;

  }
  modelprobs = REAL(Rmodelprobs);
  for(m=0; m < M; m++) {
    model = INTEGER(VECTOR_ELT(modelspace,m));
    for (j = 0; j < INTEGER(modeldim)[m]; j ++) {
      REAL(Rprobs)[model[j]] += modelprobs[m];
    }
  }
  
  for (j = 0; j < p; j ++) {
    REAL(Rprobs)[j] += (1.0 - REAL(samplemargs)[j])* eta/NC;
    if (REAL(Rprobs)[j] > 1.0) REAL(Rprobs)[j] = 1.0;
//    Rprintf("j = %d sample pip %lf corrected pip = %lf \n", j, REAL(samplemargs)[j], REAL(Rprobs)[j]);
    }
  UNPROTECT(1);
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

void model_to_vec(int *model, int p, SEXP Rmodel) {
  int j;
  memset(model, 0, p * sizeof(int));
  
  for (j = 0; j < LENGTH(Rmodel); j++) {
   model[INTEGER(Rmodel)[j]] = 1;
  }
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
      prior = exp(dpois(modeldim - 1, hyper[0], 1) - ppois(hyper[1], hyper[0], 1, 1) - lchoose((double) p-1, (double) modeldim - 1)); 
    }
  else {prior = 0.0;}

  return(prior);
}

double trunc_power_prior(int modeldim, int p, double *hyper) {
  /* modeldim and p include the intercept so subtract 1 from each */

  double prior = 0.0;
  if ((double) (modeldim -1) <= hyper[1]) {
    prior = exp(-((double) modeldim - 1.0)*((double) hyper[0])*log((double) hyper[1]+1) -
                  (log1mexp((double)(hyper[1]+1)*(double) hyper[0]*log((double) (hyper[1]+1))) - log1mexp((double) hyper[0]*log( (double) (hyper[1]+1)))) -
                  lchoose((double) p-1, (double) modeldim - 1));
    }
  else {prior = 0.0;}

  return(prior);
}

