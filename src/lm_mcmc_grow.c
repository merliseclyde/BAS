// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include "bas.h"


// [[register]]
SEXP mcmc_grow(SEXP Y, SEXP X, SEXP Rweights, SEXP Rprobinit, SEXP RnModels,
              SEXP incint, SEXP Ralpha, SEXP method, SEXP modelprior, SEXP Rupdate,
              SEXP Rbestmodel, SEXP plocal, SEXP BURNIN_Iterations,
              SEXP MCMC_Iterations, SEXP LAMBDA, SEXP DELTA,
              SEXP Rthin, SEXP Rparents, SEXP Rpivot, SEXP Rtol, SEXP Rexpand)
{

  int nModels0 = INTEGER(RnModels)[0];  // initial guess on number of models to return
  int nModels = nModels0;
  int *counts;
  

  
  double expand = REAL(Rexpand)[0]; // increase to grow vectors 
  

	// allocate return objects
	int nProtected = 0;
	
	SEXP ANS = PROTECT(allocVector(VECSXP, 16)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 16)); ++nProtected;

	SEXP Rprobs = duplicate(Rprobinit); 
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));
		
  SEXP modelspace = allocVector(VECSXP, nModels); 
  SET_VECTOR_ELT(ANS, 1, modelspace);
  SET_STRING_ELT(ANS_names, 1, mkChar("which"));
  
  SEXP Rlogmarg = allocVector(REALSXP, nModels); 
  memset(REAL(Rlogmarg), 0, nModels * sizeof(double));
  SET_VECTOR_ELT(ANS, 2, Rlogmarg);
  SET_STRING_ELT(ANS_names, 2, mkChar("logmarg"));

  SEXP modelprobs = allocVector(REALSXP, nModels);  
  memset(REAL(modelprobs), 0, nModels * sizeof(double));
  SET_VECTOR_ELT(ANS, 3, modelprobs);
  SET_STRING_ELT(ANS_names, 3, mkChar("postprobs"));
  
  SEXP priorprobs = allocVector(REALSXP, nModels); 
  memset(REAL(priorprobs), 0, nModels * sizeof(double));
  SET_VECTOR_ELT(ANS, 4, priorprobs);
  SET_STRING_ELT(ANS_names, 4, mkChar("priorprobs"));

  SEXP sampleprobs = allocVector(REALSXP, nModels); 
  memset(REAL(sampleprobs), 0, nModels * sizeof(double));
  SET_VECTOR_ELT(ANS, 5, sampleprobs);
  SET_STRING_ELT(ANS_names, 5, mkChar("sampleprobs"));

  SEXP mse = allocVector(REALSXP, nModels); 
  memset(REAL(mse), 0, nModels * sizeof(double));
  SET_VECTOR_ELT(ANS, 6, mse);
  SET_STRING_ELT(ANS_names, 6, mkChar("mse"));

  SEXP beta = allocVector(VECSXP, nModels); 
  SET_VECTOR_ELT(ANS, 7, beta);
  SET_STRING_ELT(ANS_names, 7, mkChar("mle"));
  
  SEXP se = allocVector(VECSXP, nModels);
  SET_VECTOR_ELT(ANS, 8, se);
  SET_STRING_ELT(ANS_names, 8, mkChar("mle.se"));

  SEXP shrinkage = allocVector(REALSXP, nModels); 
  memset(REAL(shrinkage), 0, nModels * sizeof(double));
  SET_VECTOR_ELT(ANS, 9, shrinkage);
  SET_STRING_ELT(ANS_names, 9, mkChar("shrinkage"));

  SEXP modeldim =  allocVector(INTSXP, nModels); 
  memset(INTEGER(modeldim), 0, nModels * sizeof(int));
  SET_VECTOR_ELT(ANS, 10, modeldim);
  SET_STRING_ELT(ANS_names, 10, mkChar("size"));
  
  SEXP R2 = allocVector(REALSXP, nModels); 
  memset(REAL(R2), 0, nModels * sizeof(double));
  SET_VECTOR_ELT(ANS, 11, R2);
  SET_STRING_ELT(ANS_names, 11, mkChar("R2"));
  
  SEXP rank = allocVector(INTSXP, nModels); 
  memset(INTEGER(rank), 0, nModels * sizeof(int));
  SET_VECTOR_ELT(ANS, 12, rank);
  SET_STRING_ELT(ANS_names, 12, mkChar("rank"));

  SEXP Rcounts =  allocVector(INTSXP, nModels); 
  counts = INTEGER(Rcounts);
  memset(counts, 0, nModels * sizeof(int));
  SET_VECTOR_ELT(ANS, 13, Rcounts);
  SET_STRING_ELT(ANS_names, 13, mkChar("freq"));
  
  SEXP MCMCprobs= duplicate(Rprobinit);
  memset(REAL(MCMCprobs), 0, LENGTH(MCMCprobs) * sizeof(double));
  SET_VECTOR_ELT(ANS, 14, MCMCprobs);
  SET_STRING_ELT(ANS_names, 14, mkChar("probne0.MCMC"));
  
  SEXP NumUnique = allocVector(INTSXP, 1); 
  SET_VECTOR_ELT(ANS, 15, NumUnique);
  SET_STRING_ELT(ANS_names, 15, mkChar("n.Unique"));
  
  setAttrib(ANS, R_NamesSymbol, ANS_names);
  
  

// allocate space for working vectors  

	SEXP RXwork = PROTECT(duplicate(X)); nProtected++;
	SEXP RYwork = PROTECT(duplicate(Y)); nProtected++;
	
	int pivot = LOGICAL(Rpivot)[0];
	double tol = REAL(Rtol)[0];

	
	double *Xwork, *Ywork,*wts, *probs, shrinkage_m,
		mse_m, MH=0.0, prior_m=1.0,
		R2_m, RSquareFull, logmarg_m, postold, postnew;
	int i, m, n, pmodel_old, *model_m, *bestmodel, rank_m;
	int mcurrent, n_sure;


	//get dimsensions of all variables
	int nobs = LENGTH(Y);
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
  //	int k = LENGTH(modelprobs);
	//	double lambda=REAL(LAMBDA)[0];
	//	double delta = REAL(DELTA)[0];
	double alpha = REAL(Ralpha)[0];
	int thin = INTEGER(Rthin)[0];


	Ywork = REAL(RYwork);
	Xwork = REAL(RXwork);
	wts = REAL(Rweights);


	double *XtXwork, *XtYwork,*XtX, *XtY, yty,SSY;
	PrecomputeData(Xwork, Ywork, wts, &XtXwork, &XtYwork, &XtX, &XtY, &yty, &SSY, p, nobs);


	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables.
	probs =  REAL(Rprobs);
	n = sortvars(vars, probs, p);
	
	for (i=n; i<p; i++) REAL(MCMCprobs)[vars[i].index] = probs[vars[i].index];
	for (i=0; i<n; i++) REAL(MCMCprobs)[vars[i].index] = 0.0;
	
	int noInclusionIs1 = no_prior_inclusion_is_1(p, probs);


	//  allocate working model and fill in the sure things
	int *model = ivecalloc(p);
	memset(model, 0, p * sizeof(int));

	for (i = n, n_sure = 0; i < p; i++)  {
		model[vars[i].index] = (int) vars[i].prob;
		if (model[vars[i].index] == 1) ++n_sure;
	}

	SEXP Rse_m = NULL, Rcoef_m = NULL, Rmodel_m=NULL;
	RSquareFull = CalculateRSquareFull(XtY, XtX, XtXwork, XtYwork, Rcoef_m, Rse_m,
                                    p, nobs, yty, SSY);

	GetRNGstate();

	NODEPTR tree, branch;
	tree = make_node(-1.0);
	//  Rprintf("For m=0, Initialize Tree with initial Model\n");

	m = 0;
	bestmodel = INTEGER(Rbestmodel);
	INTEGER(modeldim)[m] = n_sure;

	// Rprintf("Create Tree\n");
	branch = tree;
	CreateTree(branch, vars, bestmodel, model, n, m, modeldim, Rparents);

	int pmodel = INTEGER(modeldim)[m];
	PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
	memset(INTEGER(Rmodel_m), 0, pmodel * sizeof(int));
	PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
	PROTECT(Rse_m = NEW_NUMERIC(pmodel));

	model_m = GetModel_m(Rmodel_m, model, p);
	//evaluate logmarg_m and shrinkage

	R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel,
                 p, nobs, m, &mse_m, &rank_m, pivot, tol);
	INTEGER(rank)[0] = rank_m;

	gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmarg_m,
               &shrinkage_m);


	prior_m  = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
	if (prior_m == 0.0)  error("initial model has 0 prior probabilty\n");
//	SetModel2(logmarg_m, shrinkage_m, prior_m, sampleprobs, Rlogmarg, shrinkage, priorprobs, m);
  SetModel_lm(logmarg_m, shrinkage_m, prior_m, sampleprobs, Rlogmarg, shrinkage, priorprobs,
              Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2, m);

	int nUnique=1, newmodel=0, nsamples=0;
	double *real_model = vecalloc(n);
	int *modelold = ivecalloc(p);
	int old_loc = 0;
	int new_loc;
	pmodel_old = pmodel;
	
	INTEGER(Rcounts)[0] = 1;
	postold =  REAL(Rlogmarg)[m] + log(REAL(priorprobs)[m]);
	memcpy(modelold, model, sizeof(int)*p);
	m = 0;
	int *varin= ivecalloc(p);
	int *varout= ivecalloc(p);
	double problocal = REAL(plocal)[0];
	
	
	while (m < (INTEGER(MCMC_Iterations)[0] + INTEGER(BURNIN_Iterations)[0])) {

	  memcpy(model, modelold, sizeof(int)*p);
		pmodel =  n_sure;

		MH = GetNextModelCandidate(pmodel_old, n, n_sure, model, vars, problocal,
                             varin, varout, Rparents);

		branch = tree;
		newmodel= 0;
		for (i = 0; i< n; i++) {
		  int bit =  model[vars[i].index];
		  if (bit == 1) {
		    if (branch->one != NULL) branch = branch->one;
		    else newmodel = 1;
		  } else {
		    if (branch->zero != NULL)  branch = branch->zero;
		    else newmodel = 1;
		  }
		  pmodel  += bit;
		}

		if (pmodel  == n_sure || pmodel == n + n_sure) {
		  MH = 1.0/(1.0 - problocal);
		}

		if (newmodel == 1) {
		  prior_m = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
		  if (prior_m == 0.0) {
		    MH *= 0.0;
		  }
		  else {
		    new_loc = nUnique;
		    PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
		    PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
		    PROTECT(Rse_m = NEW_NUMERIC(pmodel));
		    model_m = GetModel_m(Rmodel_m, model, p);

		    R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel, p, nobs, m, &mse_m,
                      &rank_m, pivot, tol);
		    gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmarg_m, &shrinkage_m);

		    postnew = logmarg_m + log(prior_m);
		    MH *= exp(postnew - postold);
		  }}
		else {
		  new_loc = branch->where;
		  postnew =  REAL(Rlogmarg)[new_loc] +
		             log(REAL(priorprobs)[new_loc]);
		  MH *=  exp(postnew - postold);
		}

//    Rprintf("MH new %lf old %lf\n", postnew, postold);
		if (unif_rand() < MH) {
		  if (newmodel == 1) {
		    if (m % thin == 0 )  {

		    new_loc = nUnique;
		    insert_model_tree(tree, vars, n, model, nUnique);
		    INTEGER(modeldim)[nUnique] = pmodel;
		    INTEGER(rank)[nUnique] = rank_m;
		    INTEGER(Rcounts)[nUnique] = 0;  // initialize 
		    
		    //record model data
		    SetModel_lm(logmarg_m, shrinkage_m, prior_m, sampleprobs, Rlogmarg, shrinkage, priorprobs,
                    Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2,nUnique);

		    ++nUnique;
		    }
		    else UNPROTECT(3);
		  }

		  old_loc = new_loc;
		  postold = postnew;
		  pmodel_old = pmodel;
		  memcpy(modelold, model, sizeof(int)*p);

		} else  {
		  if (newmodel == 1 && prior_m > 0) UNPROTECT(3);
		}

		if ( (m % thin) == 0) {

		  INTEGER(Rcounts)[old_loc] += 1;

		  for (i = 0; i < n; i++) {
		    // store in opposite order so nth variable is first
		    real_model[n-1-i] = (double) modelold[vars[i].index];
		    REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
		  }
		  nsamples++;
		}
		if (nUnique >= nModels && m < (INTEGER(MCMC_Iterations)[0] + INTEGER(BURNIN_Iterations)[0] + 1)){
		  // expand nModels and grow result vectors
		  nModels = (int) (expand*nModels); //add checks to ensure it is not above max int

//	  Rprintf("Grow vectors:  Number of unique models %d; nModels is now %d\n", nUnique, nModels); // Need to use growable vector here

		  modelspace = resizeVector(modelspace, nModels);
		  SET_VECTOR_ELT(ANS, 1, modelspace);

		  Rlogmarg = resizeVector(Rlogmarg, nModels);
		  SET_VECTOR_ELT(ANS, 2, Rlogmarg);

		  modelprobs = resizeVector(modelprobs, nModels);
		  SET_VECTOR_ELT(ANS, 3, modelprobs);

      priorprobs = 	resizeVector(priorprobs, nModels);
      SET_VECTOR_ELT(ANS, 4, priorprobs);
      
      sampleprobs = resizeVector(sampleprobs, nModels);
		  SET_VECTOR_ELT(ANS, 5, sampleprobs);
		  
		  mse = resizeVector(mse, nModels);
		  SET_VECTOR_ELT(ANS, 6, mse);
		  
		  beta = resizeVector(beta, nModels);
		  SET_VECTOR_ELT(ANS, 7, beta);
		  
		  se = resizeVector(se, nModels);
		  SET_VECTOR_ELT(ANS, 8, se);
		  
		  shrinkage = resizeVector(shrinkage, nModels);
		  SET_VECTOR_ELT(ANS, 9, shrinkage);
		  
		  modeldim = resizeVector(modeldim, nModels);
		  SET_VECTOR_ELT(ANS, 10, modeldim);
		  
		  R2 = resizeVector(R2, nModels);
		  SET_VECTOR_ELT(ANS, 11, R2);
		  
		  rank = resizeVector(rank, nModels);
		  SET_VECTOR_ELT(ANS, 12, rank);
		  
		  Rcounts = resizeVector(Rcounts, nModels);
		  SET_VECTOR_ELT(ANS, 13, Rcounts);
		}
		m++;
	}

	// Now wrap up

	// Compute MCMC inclusion probabilities
	for (i = 0; i < n; i++) {
		REAL(MCMCprobs)[vars[i].index] /= (double) nsamples;
	}

	// Compute marginal probabilities
	mcurrent = nUnique;
	compute_modelprobs(modelprobs, Rlogmarg, priorprobs, mcurrent);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);

	INTEGER(NumUnique)[0] = nUnique;
	SET_VECTOR_ELT(ANS, 0, Rprobs);

//	Rprintf("Decreasing nModels %d to number of unique models accepted %d \n", nModels, nUnique);
	if (nUnique < nModels) {
	  SET_VECTOR_ELT(ANS, 1, resizeVector(modelspace, nUnique));
	  SET_VECTOR_ELT(ANS, 2, resizeVector(Rlogmarg, nUnique));
	  SET_VECTOR_ELT(ANS, 3, resizeVector(modelprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 4, resizeVector(priorprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 5, resizeVector(sampleprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 6, resizeVector(mse, nUnique));
	  SET_VECTOR_ELT(ANS, 7, resizeVector(beta, nUnique));
	  SET_VECTOR_ELT(ANS, 8, resizeVector(se, nUnique));
	  SET_VECTOR_ELT(ANS, 9, resizeVector(shrinkage, nUnique));
	  SET_VECTOR_ELT(ANS, 10, resizeVector(modeldim, nUnique));
	  SET_VECTOR_ELT(ANS, 11, resizeVector(R2, nUnique));
	  SET_VECTOR_ELT(ANS, 12, resizeVector(rank, nUnique));
	  SET_VECTOR_ELT(ANS, 13, resizeVector(Rcounts, nUnique));
	}	  

	PutRNGstate();
  UNPROTECT(nProtected);
    //	Rprintf("Return\n");
	return(ANS);
}

