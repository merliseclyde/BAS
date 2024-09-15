#include "bas.h"



SEXP glm_sampleworep(SEXP Y, SEXP X, SEXP Roffset, SEXP Rweights,
		     SEXP Rprobinit, SEXP Rmodeldim,
		     SEXP modelprior, SEXP betaprior,SEXP Rbestmodel,  SEXP plocal,
		     SEXP family, SEXP Rcontrol,
		     SEXP Rupdate, SEXP Rlaplace, SEXP Rparents) {
	int nProtected = 0;

	int nModels=LENGTH(Rmodeldim);

	//  Rprintf("Allocating Space for %d Models\n", nModels) ;


	SEXP ANS = PROTECT(allocVector(VECSXP, 14)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 14)); ++nProtected;
	SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
	SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP modeldim =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
	SEXP beta = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP se = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP deviance = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP priorprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP logmarg = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP sampleprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP Q = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP Rintercept = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;

	PROTECT_INDEX R2_idx;
	PROTECT_WITH_INDEX(R2, &R2_idx);
	PROTECT_INDEX shrinkage_idx;
	PROTECT_WITH_INDEX(shrinkage, &shrinkage_idx);
	PROTECT_INDEX modelspace_idx;
	PROTECT_WITH_INDEX(modelspace, &modelspace_idx);
	PROTECT_INDEX modeldim_idx;
	PROTECT_WITH_INDEX(modeldim, &modeldim_idx);
//	PROTECT_INDEX rank_idx;
//	PROTECT_WITH_INDEX(rank, &rank_idx);
	PROTECT_INDEX beta_idx;
	PROTECT_WITH_INDEX(beta, &beta_idx);
	PROTECT_INDEX se_idx;
	PROTECT_WITH_INDEX(se, &se_idx);
	PROTECT_INDEX deviance_idx;
	PROTECT_WITH_INDEX(deviance, &deviance_idx);
	PROTECT_INDEX modelprobs_idx;
	PROTECT_WITH_INDEX(modelprobs, &modelprobs_idx);
	PROTECT_INDEX priorprobs_idx;
	PROTECT_WITH_INDEX(priorprobs, &priorprobs_idx);
	PROTECT_INDEX logmarg_idx;
	PROTECT_WITH_INDEX(logmarg, &logmarg_idx);
	PROTECT_INDEX sampleprobs_idx;
	PROTECT_WITH_INDEX(sampleprobs, &sampleprobs_idx);
	PROTECT_INDEX Q_idx;
	PROTECT_WITH_INDEX(R2, &Q_idx);
	PROTECT_INDEX Rintercept_idx;
	PROTECT_WITH_INDEX(Rintercept, &Rintercept_idx);
	
	
	double *probs,logmargy, shrinkage_m;
	int i;

	glmstptr *glmfamily;
	glmfamily = make_glmfamily_structure(family);

	betapriorptr *betapriorfamily;
	betapriorfamily = make_betaprior_structure(betaprior, family);


	//get dimsensions of all variables
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int k = LENGTH(modelprobs);

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

	SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
	SetModel1(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q,Rintercept, m);

	UNPROTECT(2);

	int *modelwork= ivecalloc(p);

	// sample models
	for (m = 1;  m < k  && lessThanOne(pigamma[0]); m++) {
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

		SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
		SetModel1(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2,Q,Rintercept, m);
		UNPROTECT(2);

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
				if (update_probs(probs, vars, mcurrent, k, p) == 1) {
				  //					Rprintf("Updating Model Tree %d \n", m);
					update_tree(modelspace, tree, modeldim, vars, k,p,n,mcurrent, modelwork);
				}
			}
		}
	}

	if (m < k) {
	  // resize if constraints have reduced the number of models
	  k = m;
	  
	  REPROTECT(logmarg= Rf_lengthgets(logmarg, m), logmarg_idx);
	  REPROTECT(modelprobs= Rf_lengthgets(modelprobs, m), modelprobs_idx);
	  REPROTECT(priorprobs= Rf_lengthgets(priorprobs, m), priorprobs_idx);
	  REPROTECT(sampleprobs= Rf_lengthgets(sampleprobs, m), sampleprobs_idx);
	  REPROTECT(deviance = Rf_lengthgets(deviance, m), deviance_idx);
	  REPROTECT(shrinkage = Rf_lengthgets(shrinkage, m), shrinkage_idx);
	  REPROTECT(modeldim= Rf_lengthgets(modeldim, m), modeldim_idx);
	  REPROTECT(R2= Rf_lengthgets(R2, m), R2_idx);
	  REPROTECT(se= Rf_lengthgets(se, m), se_idx);
//	  REPROTECT(rank = Rf_lengthgets(rank, m), rank_idx);
	  REPROTECT(modelspace = Rf_lengthgets(modelspace, m), modelspace_idx);
	  REPROTECT(beta = Rf_lengthgets(beta, m), beta_idx);
	  REPROTECT(se= Rf_lengthgets(se, m), se_idx);
	  REPROTECT(Q= Rf_lengthgets(Q, m), Q_idx);
	  REPROTECT(Rintercept= Rf_lengthgets(Rintercept, m), Rintercept_idx);

	}

	compute_modelprobs(modelprobs, logmarg, priorprobs,k);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, k, p);

	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));

	SET_VECTOR_ELT(ANS, 1, modelspace);
	SET_STRING_ELT(ANS_names, 1, mkChar("which"));

	SET_VECTOR_ELT(ANS, 2, logmarg);
	SET_STRING_ELT(ANS_names, 2, mkChar("logmarg"));

	SET_VECTOR_ELT(ANS, 3, modelprobs);
	SET_STRING_ELT(ANS_names, 3, mkChar("postprobs"));

	SET_VECTOR_ELT(ANS, 4, priorprobs);
	SET_STRING_ELT(ANS_names, 4, mkChar("priorprobs"));

	SET_VECTOR_ELT(ANS, 5,sampleprobs);
	SET_STRING_ELT(ANS_names, 5, mkChar("sampleprobs"));

	SET_VECTOR_ELT(ANS, 6, deviance);
	SET_STRING_ELT(ANS_names, 6, mkChar("deviance"));

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

	SET_VECTOR_ELT(ANS, 12, Q);
	SET_STRING_ELT(ANS_names, 12, mkChar("Q"));

	SET_VECTOR_ELT(ANS, 13, Rintercept);
	SET_STRING_ELT(ANS_names, 13, mkChar("intercept"));

	setAttrib(ANS, R_NamesSymbol, ANS_names);
	PutRNGstate();

	UNPROTECT(nProtected);
	return(ANS);
}
