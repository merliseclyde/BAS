// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
/* version  5/20/2005 */
/* Rsample.c program for sampling without replacement in R  MC 11/2002 */
/* based on sim.c: program for running simulations with random and
deterministic sampling. ML 6/97. */
/*  top-k.c: Michael Littman, Sun Dec 15 19:29:05 EST 1996
*   Version 4.  Assume entries are positive and sorted (big to small).
*  Given a set of n integers, list the k subsets that have the
*  highest sums (in order).
*
* Michael Littman, Tue Jun  3 11:38:08 EDT 1997
*  Modifying to run more standalone.  In particular, does the logit
*  calculations and sorting itself instead of depending on S to do it.
* Merlise Clyde, February 2003,  modified to be called from R
* reworked memory management and tree structures for larger problems
*/

/* Includes. */
#include "bas.h"

// extern inline int lessThanOne(double a);

// [[register]]
extern SEXP sampleworep_new(SEXP Y, SEXP X, SEXP Rweights, SEXP Rprobinit,
                            SEXP Rmodeldim, SEXP incint, SEXP Ralpha,
                            SEXP method, SEXP modelprior, SEXP Rupdate,
                            SEXP Rbestmodel, SEXP plocal, SEXP Rparents, 
                            SEXP Rpivot, SEXP Rtol) {
	int nProtected = 0;
	SEXP RXwork = PROTECT(duplicate(X)); nProtected++;
	SEXP RYwork = PROTECT(duplicate(Y)); nProtected++;
	int nModels=LENGTH(Rmodeldim);
	int pivot = LOGICAL(Rpivot)[0];
	double tol = REAL(Rtol)[0];
	
	
	int nModels0 = LENGTH(Rmodeldim);;  // initial guess on number of models to return
	int nUnique = nModels0;

	SEXP ANS = PROTECT(allocVector(VECSXP, 14)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 14)); ++nProtected;
	
	SEXP Rprobs = duplicate(Rprobinit); 
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));
	
	SEXP modelspace = allocVector(VECSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 1, modelspace);
	SET_STRING_ELT(ANS_names, 1, mkChar("which"));
	
	SEXP Rlogmarg = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 2, Rlogmarg);
	SET_STRING_ELT(ANS_names, 2, mkChar("logmarg"));
	
	SEXP modelprobs = allocVector(REALSXP, nModels0);  
	SET_VECTOR_ELT(ANS, 3, modelprobs);
	SET_STRING_ELT(ANS_names, 3, mkChar("postprobs"));
	
	SEXP priorprobs = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 4, priorprobs);
	SET_STRING_ELT(ANS_names, 4, mkChar("priorprobs"));
	
	SEXP sampleprobs = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 5, sampleprobs);
	SET_STRING_ELT(ANS_names, 5, mkChar("sampleprobs"));
	
	SEXP mse = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 6, mse);
	SET_STRING_ELT(ANS_names, 6, mkChar("mse"));
	
	SEXP beta = allocVector(VECSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 7, beta);
	SET_STRING_ELT(ANS_names, 7, mkChar("mle"));
	
	SEXP se = allocVector(VECSXP, nModels0);
	SET_VECTOR_ELT(ANS, 8, se);
	SET_STRING_ELT(ANS_names, 8, mkChar("mle.se"));
	
	SEXP shrinkage = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 9, shrinkage);
	SET_STRING_ELT(ANS_names, 9, mkChar("shrinkage"));
	
	SEXP modeldim =  allocVector(INTSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 10, modeldim);
	SET_STRING_ELT(ANS_names, 10, mkChar("size"));
	
	SEXP R2 = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 11, R2);
	SET_STRING_ELT(ANS_names, 11, mkChar("R2"));
	
	SEXP rank = allocVector(INTSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 12, rank);
	SET_STRING_ELT(ANS_names, 12, mkChar("rank"));
	
	
	SEXP NumUnique = allocVector(INTSXP, 1); 
	SET_VECTOR_ELT(ANS, 13, NumUnique);
	SET_STRING_ELT(ANS_names, 13, mkChar("n.Unique"));
	
	setAttrib(ANS, R_NamesSymbol, ANS_names);
	
	memset(REAL(Rlogmarg), 0.0, sizeof(double) *nModels);
	memset(REAL(modelprobs), 0.0, sizeof(double) *nModels);
	memset(REAL(priorprobs), 0.0, sizeof(double) *nModels);
	memset(REAL(shrinkage), 0.0, sizeof(double) *nModels);
	memset(REAL(R2), 0.0, sizeof(double) *nModels);
	memset(REAL(Rlogmarg), 0.0, sizeof(double) *nModels);
	memset(REAL(sampleprobs), 0.0, sizeof(double) *nModels);
	memset(INTEGER(modeldim), 0, sizeof(int) *nModels);
	memset(REAL(mse), 0.0, sizeof(double) *nModels);
	memset(INTEGER(rank), 0, sizeof(int) *nModels);	

	//  Rprintf("Allocating Space for %d Models\n", nModels) ;
	/* Allocate Space for Results  R < 4.5.2
	SEXP ANS = PROTECT(allocVector(VECSXP, 13)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 13)); ++nProtected;
	SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
	

	  SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	  SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	  SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	  SEXP modeldim =  PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	  SEXP rank =  PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	  SEXP beta = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	  SEXP se = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	  SEXP mse = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	  SEXP modelprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	  SEXP priorprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	  SEXP Rlogmarg = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	  SEXP sampleprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
*/

	double *Xwork, *Ywork, *wts, *probs, shrinkage_m, mse_m, R2_m, RSquareFull, Rbestmarg, logmarg_m;
	int i, *model_m, *bestmodel, rank_m;

	//get dimsensions of all variables
	int nobs = LENGTH(Y);
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int k = LENGTH(modelprobs);
	double alpha = REAL(Ralpha)[0];
	SEXP Rbestmodel_new = PROTECT(duplicate(Rbestmodel)); nProtected++;

	int update = INTEGER(Rupdate)[0];
	double eps = DBL_EPSILON;
	double problocal = REAL(plocal)[0];

  // memset(INTEGER(modeldim), 1, k*sizeof(int));
	Ywork = REAL(RYwork);
	Xwork = REAL(RXwork);
	wts = REAL(Rweights);

	double *XtXwork, *XtYwork,*XtX, *XtY, yty=0.0,SSY=0.0;
	PrecomputeData(Xwork, Ywork, wts, &XtXwork, &XtYwork, &XtX, &XtY, &yty, &SSY, p, nobs);

	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables.
	probs =  REAL(Rprobs);
	int n = sortvars(vars, probs, p);
	int noInclusionIs1 = no_prior_inclusion_is_1(p, probs);

	SEXP  Rse_m = NULL, Rcoef_m = NULL, Rmodel_m = NULL;
	RSquareFull = CalculateRSquareFull(XtY, XtX, XtXwork, XtYwork, Rcoef_m, Rse_m, p, nobs, yty, SSY);

	int *model = ivecalloc(p);
	memset(model, 0, p * sizeof(int));

	/* fill in the sure things */
	for (i = n; i < p; i++)  {
		model[vars[i].index] = (int) vars[i].prob;
	}

	GetRNGstate();

	NODEPTR tree, branch;
	tree = make_node(vars[0].prob);
	//	Rprintf("For m=0, Initialize Tree with initial Model\n");

	int m = 0;
	bestmodel = INTEGER(Rbestmodel_new);
	REAL(Rlogmarg)[m] = 0.0;
	INTEGER(modeldim)[m] = 0;

	for (i = n; i < p; i++)  {
		model[vars[i].index] = bestmodel[vars[i].index];
		INTEGER(modeldim)[m]  +=  bestmodel[vars[i].index];
	}
	double *pigamma = vecalloc(p);
	memset(pigamma, 0.0 ,p*sizeof(double)); 
	branch = tree;
	CreateTree_with_pigamma(branch, vars, bestmodel, model, n, m,
                         modeldim, pigamma, Rparents);

	branch=tree;
	Substract_visited_probability_mass(branch, vars, model, n, m, pigamma,eps);

	int pmodel = INTEGER(modeldim)[m];


	PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
	memset(INTEGER(Rmodel_m), 0, pmodel * sizeof(int));
	PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
	PROTECT(Rse_m = NEW_NUMERIC(pmodel));

  model_m = GetModel_m(Rmodel_m, model, p);

	R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY,
                 pmodel, p, nobs, m, &mse_m, &rank_m, pivot, tol);
	INTEGER(rank)[m] = rank_m;

	gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmarg_m, &shrinkage_m);
//	Rprintf("rank %d dim %d\n", rank_m, pmodel);
//	gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmarg_m, &shrinkage_m);

//  check should this depend on rank or pmodel?
	double prior_m  = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);



//	SetModel2(logmarg_m, shrinkage_m, prior_m, sampleprobs, Rlogmarg, shrinkage, priorprobs, m);
	SetModel_lm(logmarg_m, shrinkage_m, prior_m, sampleprobs, Rlogmarg, shrinkage, priorprobs,
             Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2, m);
	//Rprintf("model %d max logmarg %lf\n", m, REAL(Rlogmarg)[m]);

	Rbestmarg = REAL(Rlogmarg)[m];
//  double *parents = REAL(Rparents);
	int *modelwork= ivecalloc(p);
	/*
  for (j =0; j < p; j++) Rprintf("%d ", vars[j].index);
    Rprintf("\n");
	for (i=0; i < p; i++) {
	  Rprintf("%d ", vars[i].index);
	  for (j = 0; j < p; j++) {
	    Rprintf("%lf ", parents[vars[i].index + p*vars[j].index]);
	  }
	  Rprintf("\n");
	}  */

	// Sample models
	//	for (m = 1;  m < k && pigamma[0] < 1.0; m++) {
	for (m = 1;  m < k && lessThanOne(pigamma[0]); m++) {
	//  Rprintf("model %d, starting pigamma = %lf\n", m, pigamma[0]);
	  INTEGER(modeldim)[m] = 0;
		for (i = n; i < p; i++)  {
			INTEGER(modeldim)[m]  +=  model[vars[i].index];
		}

		branch = tree;
		GetNextModel_swop(branch, vars, model, n, m, pigamma, problocal,
                    modeldim, bestmodel, Rparents);

		/* Now subtract off the visited probability mass. */
		branch=tree;
		Substract_visited_probability_mass(branch, vars, model, n, m, pigamma,eps);

		/* Now get model specific calculations */
		pmodel = INTEGER(modeldim)[m];
		PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
		memset(INTEGER(Rmodel_m), 0, pmodel * sizeof(int));
		PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
		PROTECT(Rse_m = NEW_NUMERIC(pmodel));
		model_m = GetModel_m(Rmodel_m, model, p);

		R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY,
                  pmodel, p, nobs, m, &mse_m, &rank_m, pivot, tol);
		INTEGER(rank)[m] = rank_m;
		// initialize
		logmarg_m= 0.0;
		shrinkage_m = 1.0;
		gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmarg_m, &shrinkage_m);
//    Rprintf("rank %d dim %d\n", rank_m, pmodel);
//		gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmarg_m, &shrinkage_m);

		prior_m = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
//		SetModel2(logmarg_m, shrinkage_m, prior_m, sampleprobs, Rlogmarg, shrinkage, priorprobs, m);
		SetModel_lm(logmarg_m, shrinkage_m, prior_m, sampleprobs, Rlogmarg, shrinkage, priorprobs,
                Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2,m);
	  

		REAL(sampleprobs)[m] = pigamma[0];

		//update best model
		if (REAL(Rlogmarg)[m] > Rbestmarg) {
			for (i=0; i < p; i++) {
				bestmodel[i] = model[i];
			}
			Rbestmarg = REAL(Rlogmarg)[m];
		}

		//update marginal inclusion probs
		if (m > 1) {
			double mod;
			double rem = modf((double) m/(double) update, &mod);
			if (rem  == 0.0) {
				int mcurrent = m;
				compute_modelprobs(modelprobs, Rlogmarg, priorprobs,mcurrent);
				compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);
				if (update_probs(probs, vars, mcurrent, k, p) == 1) {
				  //					Rprintf("Updating Model Tree %d \n", m);
					update_tree(modelspace, tree, modeldim, vars, k,p,n,mcurrent, modelwork);
				}
			}
		}
	}
  nUnique = k;
	if (m < nUnique) { 
	  // Rprintf("allocated %d models but only %d sampled; resizing\n", nUnique, m); 
	  // resize
	  // consider using force.heredity
	  nUnique = m;
	  
	  compute_modelprobs(modelprobs, Rlogmarg, priorprobs,nUnique);
	  compute_margprobs(modelspace, modeldim, modelprobs, probs, nUnique, p);
	  SET_VECTOR_ELT(ANS, 1, xlengthgets(modelspace, nUnique));
	  SET_VECTOR_ELT(ANS, 2, xlengthgets(Rlogmarg, nUnique));
	  SET_VECTOR_ELT(ANS, 3, xlengthgets(modelprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 4, xlengthgets(priorprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 5, xlengthgets(sampleprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 6, xlengthgets(mse, nUnique));
	  SET_VECTOR_ELT(ANS, 7, xlengthgets(beta, nUnique));
	  SET_VECTOR_ELT(ANS, 8, xlengthgets(se, nUnique));
	  SET_VECTOR_ELT(ANS, 9, xlengthgets(shrinkage, nUnique));
	  SET_VECTOR_ELT(ANS, 10, xlengthgets(modeldim, nUnique));
	  SET_VECTOR_ELT(ANS, 11, xlengthgets(R2, nUnique));
	  SET_VECTOR_ELT(ANS, 12, xlengthgets(rank, nUnique));
	}
	else {
	  compute_modelprobs(modelprobs, Rlogmarg, priorprobs,nUnique);
	  compute_margprobs(modelspace, modeldim, modelprobs, probs, nUnique, p);
	}
	
	INTEGER(NumUnique)[0] = nUnique;
/*	
  if (m < k) { 
 // warning("allocated %d models but only %d sampled; using SETLENGTH to resize\n", k, m); 
 // resize
 // consider using force.heredity
    k = m;

    SETLENGTH(modelspace, m);
    SETLENGTH(Rlogmarg, m);
    SETLENGTH(modelprobs, m);
    SETLENGTH(priorprobs, m);
    SETLENGTH(sampleprobs, m);
    SETLENGTH(beta, m);
    SETLENGTH(se, m);
    SETLENGTH(mse, m);
    SETLENGTH(shrinkage, m);
    SETLENGTH(modeldim, m);
    SETLENGTH(R2, m);
    SETLENGTH(rank, m);
  }


	compute_modelprobs(modelprobs, Rlogmarg, priorprobs,k);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, k, p);

	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));

	SET_VECTOR_ELT(ANS, 1, modelspace);
	SET_STRING_ELT(ANS_names, 1, mkChar("which"));

	SET_VECTOR_ELT(ANS, 2, Rlogmarg);
	SET_STRING_ELT(ANS_names, 2, mkChar("logmarg"));

	SET_VECTOR_ELT(ANS, 3, modelprobs);
	SET_STRING_ELT(ANS_names, 3, mkChar("postprobs"));

	SET_VECTOR_ELT(ANS, 4, priorprobs);
	SET_STRING_ELT(ANS_names, 4, mkChar("priorprobs"));

	SET_VECTOR_ELT(ANS, 5,sampleprobs);
	SET_STRING_ELT(ANS_names, 5, mkChar("sampleprobs"));

	SET_VECTOR_ELT(ANS, 6, mse);
	SET_STRING_ELT(ANS_names, 6, mkChar("mse"));

	SET_VECTOR_ELT(ANS, 7, beta);
	SET_STRING_ELT(ANS_names, 7, mkChar("mle"));

	SET_VECTOR_ELT(ANS, 8, se);
	SET_STRING_ELT(ANS_names, 8, mkChar("mle.se"));

	SET_VECTOR_ELT(ANS, 9, shrinkage);
	SET_STRING_ELT(ANS_names, 9, mkChar("shrinkage"));

	SET_VECTOR_ELT(ANS, 10, modeldim);
	SET_STRING_ELT(ANS_names, 10, mkChar("size"));

	SET_VECTOR_ELT(ANS, 11, R2);
	SET_STRING_ELT(ANS_names, 11, mkChar("R2"));

	SET_VECTOR_ELT(ANS, 12, rank);
	SET_STRING_ELT(ANS_names, 12, mkChar("rank"));

	setAttrib(ANS, R_NamesSymbol, ANS_names);
 
 */
	PutRNGstate();

	UNPROTECT(nProtected);


	return(ANS);
}
