// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include "bas.h"


SEXP glm_sampleworep_grow(SEXP Y, SEXP X, SEXP Roffset, SEXP Rweights,
		     SEXP Rprobinit, SEXP RnModels,
		     SEXP modelprior, SEXP betaprior,SEXP Rbestmodel,  SEXP plocal,
		     SEXP family, SEXP Rcontrol,
		     SEXP Rupdate, SEXP Rlaplace, SEXP Rparents) {


  int nModels0 = INTEGER(RnModels)[0];  // initial guess on number of models to return
  int nUnique = nModels0;

//	Rprintf("Allocating Space for %d Models\n", nModels0) ;

	int nProtected = 0;
	
	SEXP ANS = PROTECT(allocVector(VECSXP, 15)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 15)); ++nProtected;
	
	SEXP Rprobs = duplicate(Rprobinit); 
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));
	
	SEXP modelspace = allocVector(VECSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 1, modelspace);
	SET_STRING_ELT(ANS_names, 1, mkChar("which"));
	
	SEXP logmarg = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 2, logmarg);
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
	
	SEXP deviance = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 6, deviance);
	SET_STRING_ELT(ANS_names, 6, mkChar("deviance"));
	
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
	memset(INTEGER(modeldim), 0, nModels0 * sizeof(int));
	SET_VECTOR_ELT(ANS, 10, modeldim);
	SET_STRING_ELT(ANS_names, 10, mkChar("size"));
	
	SEXP R2 = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 11, R2);
	SET_STRING_ELT(ANS_names, 11, mkChar("R2"));
	
	SEXP NumUnique = allocVector(INTSXP, 1); 
	SET_VECTOR_ELT(ANS, 12, NumUnique);
	SET_STRING_ELT(ANS_names, 12, mkChar("n.Unique"));
	
	SEXP Q = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 13, Q);
	SET_STRING_ELT(ANS_names, 13, mkChar("Q"));
	
	SEXP Rintercept = allocVector(REALSXP, nModels0); 
	SET_VECTOR_ELT(ANS, 14, Rintercept);
	SET_STRING_ELT(ANS_names, 14, mkChar("intercept"));
	
	setAttrib(ANS, R_NamesSymbol, ANS_names);
	
//	Rprintf("Start Computing\n");
	double *probs,logmargy, shrinkage_m;
	int i;

	glmstptr *glmfamily;
	glmfamily = make_glmfamily_structure(family);

	betapriorptr *betapriorfamily;
	betapriorfamily = make_betaprior_structure(betaprior, family);


	//get dimsensions of all variables
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];

	int update = INTEGER(Rupdate)[0];
	double eps = DBL_EPSILON;
	double problocal = REAL(plocal)[0];

	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables.
	probs =  REAL(Rprobs);
	int n = sortvars(vars, probs, p);
	int noInclusionIs1 = no_prior_inclusion_is_1(p, probs);

	int *model = ivecalloc(p);
	/* fill in the sure things */
	for (i = n; i < p; i++)  {
		model[vars[i].index] = (int) vars[i].prob;
	}

	GetRNGstate();

	NODEPTR tree, branch;
	tree = make_node(vars[0].prob);
//	Rprintf("For m=0, Initialize Tree with initial Model\n");

	int m = 0;
	int *bestmodel = INTEGER(Rbestmodel);
	
	INTEGER(modeldim)[m] = 0;
	for (i = n; i < p; i++)  {
		model[vars[i].index] = bestmodel[vars[i].index];
		INTEGER(modeldim)[m]  +=  bestmodel[vars[i].index];
	}

	double *pigamma = vecalloc(p);
	branch = tree;
	CreateTree_with_pigamma(branch, vars, bestmodel, model, n, m, modeldim,pigamma, Rparents);

	branch=tree;
	Substract_visited_probability_mass(branch, vars, model, n, m, pigamma,eps);

	int pmodel = INTEGER(modeldim)[m];
	SEXP Rmodel_m =	PROTECT(allocVector(INTSXP,pmodel));
	GetModel_m(Rmodel_m, model, p);
	//evaluate logmargy and shrinkage
	SEXP glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights,
					    glmfamily, Rcontrol, Rlaplace,
					    betapriorfamily));
	double prior_m  = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
	logmargy = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
	shrinkage_m = REAL(getListElement(getListElement(glm_fit, "lpy"),
					"shrinkage"))[0];
//   Rprintf("SetModel_glm for initial model\n");
//	SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
//	SetModel1(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q,Rintercept, m);
	SetModel_glm(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q,Rintercept, 
              prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
//	UNPROTECT(2);

	int *modelwork= ivecalloc(p);

//	Rprintf("sample models\n");
	for (m = 1;  m < nUnique  && lessThanOne(pigamma[0]); m++) {
	  INTEGER(modeldim)[m] = 0.0;
		for (i = n; i < p; i++)  {
			INTEGER(modeldim)[m]  +=  model[vars[i].index];
		}

		branch = tree;
		GetNextModel_swop(branch, vars, model, n, m, pigamma, problocal,
                      modeldim, bestmodel,Rparents);

		/* Now subtract off the visited probability mass. */
		branch=tree;
		Substract_visited_probability_mass(branch, vars, model, n, m, pigamma,eps);

		/* Now get model specific calculations */
		pmodel = INTEGER(modeldim)[m];
		PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
		memset(INTEGER(Rmodel_m), 0, pmodel * sizeof(int));
		GetModel_m(Rmodel_m, model, p);

		glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights,
					       glmfamily, Rcontrol, Rlaplace,
					       betapriorfamily));
		prior_m = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
		logmargy = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
		shrinkage_m = REAL(getListElement(getListElement(glm_fit, "lpy"),
					"shrinkage"))[0];

//		SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
//		SetModel1(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2,Q,Rintercept, m);
    SetModel_glm(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q,Rintercept, 
                 prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
//		UNPROTECT(2);

		REAL(sampleprobs)[m] = pigamma[0];

		//update best model

		//update marginal inclusion probs
		if (m > 1) {
			double mod;
			double rem = modf((double) m/(double) update, &mod);
			if (rem  == 0.0) {
				int mcurrent = m;
				compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
				compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);
				if (update_probs(probs, vars, mcurrent, nUnique, p) == 1) {
				  //					Rprintf("Updating Model Tree %d \n", m);
					update_tree(modelspace, tree, modeldim, vars, nUnique, p, n, mcurrent, modelwork);
				}
			}
		}
	}

	
	if (m < nUnique) {
//	  Rprintf("resize if constraints have reduced the number of models\n");
	  nUnique = m;
	  
	  SET_VECTOR_ELT(ANS, 1, resizeVector(modelspace, nUnique));
	  SET_VECTOR_ELT(ANS, 2, resizeVector(logmarg, nUnique));
	  SET_VECTOR_ELT(ANS, 3, resizeVector(modelprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 4, resizeVector(priorprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 5, resizeVector(sampleprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 6, resizeVector(deviance, nUnique));
	  SET_VECTOR_ELT(ANS, 7, resizeVector(beta, nUnique));
	  SET_VECTOR_ELT(ANS, 8, resizeVector(se, nUnique));
	  SET_VECTOR_ELT(ANS, 9, resizeVector(shrinkage, nUnique));
	  SET_VECTOR_ELT(ANS, 10, resizeVector(modeldim, nUnique));
	  SET_VECTOR_ELT(ANS, 11, resizeVector(R2, nUnique));
	  SET_VECTOR_ELT(ANS, 13, resizeVector(Q, nUnique));
	  SET_VECTOR_ELT(ANS, 14, resizeVector(Rintercept, nUnique));
//	  Rprintf("resizing to %d models\n", nUnique);
	}

	compute_modelprobs(modelprobs, logmarg, priorprobs,nUnique);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, nUnique, p);

//	Rprintf("computed model probs\n");
	INTEGER(NumUnique)[0] = nUnique;
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	
//	Rprintf("returning\n");
	
	PutRNGstate();

	UNPROTECT(nProtected);
	return(ANS);
}
