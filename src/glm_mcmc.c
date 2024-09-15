#include "bas.h"


// [[register]]
SEXP glm_mcmc(SEXP Y, SEXP X, SEXP Roffset, SEXP Rweights,
	      SEXP Rprobinit, SEXP Rmodeldim,
	      SEXP modelprior,  SEXP betaprior, SEXP Rbestmodel,  SEXP plocal,
	      SEXP BURNIN_Iterations, SEXP Rthin, 
	      SEXP family, SEXP Rcontrol, SEXP Rlaplace, SEXP Rparents
			  )
{
	int nProtected = 0;
	int nModels=LENGTH(Rmodeldim);
	SEXP ANS = PROTECT(allocVector(VECSXP, 17)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 17)); ++nProtected;
	SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
	SEXP MCMCprobs= PROTECT(duplicate(Rprobinit)); ++nProtected;
	SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP modeldim =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
	SEXP counts =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
	SEXP beta = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP se = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP deviance = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP priorprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP logmarg = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP sampleprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP Q = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP Rintercept = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;

	SEXP NumUnique = PROTECT(allocVector(INTSXP, 1)); ++nProtected;

	PROTECT_INDEX counts_idx;
	PROTECT_WITH_INDEX(counts, &counts_idx);
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
	
	
	double *probs, MH=0.0, prior_m=1.0, shrinkage_m, logmargy, postold, postnew;
	int i, m, n, pmodel_old, *bestmodel;
	int mcurrent, n_sure;

	glmstptr *glmfamily;
	glmfamily = make_glmfamily_structure(family);

	betapriorptr *betapriorfamily;
	betapriorfamily = make_betaprior_structure(betaprior, family);


	//get dimsensions of all variables
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int k = LENGTH(modelprobs);
	
	int thin = INTEGER(Rthin)[0];

	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables.
	probs =  REAL(Rprobs);
	n = sortvars(vars, probs, p);
	for (i =n; i <p; i++) REAL(MCMCprobs)[vars[i].index] = probs[vars[i].index];
	for (i =0; i <n; i++) REAL(MCMCprobs)[vars[i].index] = 0.0;
	int noInclusionIs1 = no_prior_inclusion_is_1(p, probs);

	// fill in the sure things
	int *model = ivecalloc(p);
	for (i = n, n_sure = 0; i < p; i++)  {
		model[vars[i].index] = (int) vars[i].prob;
		if (model[vars[i].index] == 1) ++n_sure;
	}

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
	SEXP Rmodel_m =	PROTECT(allocVector(INTSXP,pmodel));
	GetModel_m(Rmodel_m, model, p);
	//evaluate logmargy and shrinkage
	SEXP glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights,
					    glmfamily, Rcontrol, Rlaplace,
					    betapriorfamily));
	prior_m  = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);

	logmargy = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
	shrinkage_m = REAL(getListElement(getListElement(glm_fit, "lpy"),
					"shrinkage"))[0];
	SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
	SetModel1(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q,Rintercept, m);
	UNPROTECT(2);

	int nUnique=0, newmodel=0;
	double *real_model = vecalloc(n);
	int *modelold = ivecalloc(p);
	int old_loc = 0;
	int new_loc;
	pmodel_old = pmodel;
	nUnique=1;
	INTEGER(counts)[0] = 1;
	postold =  REAL(logmarg)[m] + log(REAL(priorprobs)[m]);
	memcpy(modelold, model, sizeof(int)*p);
	m = 0;
	int *varin= ivecalloc(p);
	int *varout= ivecalloc(p);
	double problocal = REAL(plocal)[0];
	while (nUnique < k && m < INTEGER(BURNIN_Iterations)[0]) {
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
		  new_loc = nUnique;
		  PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
		  GetModel_m(Rmodel_m, model, p);

		  glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights,
						 glmfamily, Rcontrol, Rlaplace,
						 betapriorfamily));
		  prior_m = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);

		  logmargy = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
		  shrinkage_m = REAL(getListElement(getListElement(glm_fit, "lpy"),
						  "shrinkage"))[0];

		  postnew = logmargy + log(prior_m);
		} else {
		  new_loc = branch->where;
		  postnew =  REAL(logmarg)[new_loc] + log(REAL(priorprobs)[new_loc]);
		}

		MH *= exp(postnew - postold);
		//    Rprintf("MH new %lf old %lf\n", postnew, postold);
		if (unif_rand() < MH) {
		 if (newmodel == 1)  {
			if ((m % thin) == 0 )  {
			  new_loc = nUnique;
			  insert_model_tree(tree, vars, n, model, nUnique);
			  INTEGER(modeldim)[nUnique] = pmodel;
				//Rprintf("model %d: %d variables\n", m, pmodel);
			  SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, nUnique);
			  SetModel1(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q, Rintercept, nUnique);
			  ++nUnique;
			}
			UNPROTECT(2);
		 }
			old_loc = new_loc;
			postold = postnew;
			pmodel_old = pmodel;
			memcpy(modelold, model, sizeof(int)*p);
		 } else  {
			if (newmodel == 1) UNPROTECT(2);
		}
		INTEGER(counts)[old_loc] += 1;
		for (i = 0; i < n; i++) {
			// store in opposite order so nth variable is first
			real_model[n-1-i] = (double) modelold[vars[i].index];
			REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
		}
	m++;
	}

	for (i = 0; i < n; i++) {
		REAL(MCMCprobs)[vars[i].index] /= (double) m;
	}



	// Compute marginal probabilities
	mcurrent = nUnique;
	//	Rprintf("NumUnique Models Accepted %d \n", nUnique);
	compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);

	INTEGER(NumUnique)[0] = nUnique;
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));

	if (nUnique < nModels) {
	  nModels = nUnique;
	  REPROTECT(logmarg= Rf_lengthgets(logmarg, nUnique), logmarg_idx);
	  REPROTECT(modelprobs= Rf_lengthgets(modelprobs, nUnique), modelprobs_idx);
	  REPROTECT(priorprobs= Rf_lengthgets(priorprobs, nUnique), priorprobs_idx);
	  REPROTECT(sampleprobs= Rf_lengthgets(sampleprobs, nUnique), sampleprobs_idx);
	  REPROTECT(deviance = Rf_lengthgets(deviance, nUnique), deviance_idx);
	  REPROTECT(shrinkage = Rf_lengthgets(shrinkage, nUnique), shrinkage_idx);
	  REPROTECT(modeldim= Rf_lengthgets(modeldim, nUnique), modeldim_idx);
	  REPROTECT(R2= Rf_lengthgets(R2, nUnique), R2_idx);
	  REPROTECT(se= Rf_lengthgets(se, nUnique), se_idx);
	  //	  REPROTECT(rank = Rf_lengthgets(rank, nUnique), rank_idx);
	  REPROTECT(modelspace = Rf_lengthgets(modelspace, nUnique), modelspace_idx);
	  REPROTECT(beta = Rf_lengthgets(beta, nUnique), beta_idx);
	  REPROTECT(se= Rf_lengthgets(se, nUnique), se_idx);
	  REPROTECT(Q= Rf_lengthgets(Q, nUnique), Q_idx);
	  REPROTECT(Rintercept= Rf_lengthgets(Rintercept, nUnique), Rintercept_idx);
	  REPROTECT(counts= Rf_lengthgets(counts, nUnique), counts_idx);
	  
	}

	SET_VECTOR_ELT(ANS, 1, modelspace);
	SET_STRING_ELT(ANS_names, 1, mkChar("which"));

	SET_VECTOR_ELT(ANS, 2, logmarg);
	SET_STRING_ELT(ANS_names, 2, mkChar("logmarg"));

	SET_VECTOR_ELT(ANS, 3, modelprobs);
	SET_STRING_ELT(ANS_names, 3, mkChar("postprobs"));

	SET_VECTOR_ELT(ANS, 4, priorprobs);
	SET_STRING_ELT(ANS_names, 4, mkChar("priorprobs"));

	SET_VECTOR_ELT(ANS, 5, sampleprobs);
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

	SET_VECTOR_ELT(ANS, 12, counts);
	SET_STRING_ELT(ANS_names, 12, mkChar("freq"));

	SET_VECTOR_ELT(ANS, 13, MCMCprobs);
	SET_STRING_ELT(ANS_names, 13, mkChar("probne0.MCMC"));

	SET_VECTOR_ELT(ANS, 14, NumUnique);
	SET_STRING_ELT(ANS_names, 14, mkChar("n.Unique"));

  SET_VECTOR_ELT(ANS, 15, Q);
  SET_STRING_ELT(ANS_names, 15, mkChar("Q"));

	SET_VECTOR_ELT(ANS, 16, Rintercept);
	SET_STRING_ELT(ANS_names, 16, mkChar("intercept"));

	setAttrib(ANS, R_NamesSymbol, ANS_names);

	PutRNGstate();
	UNPROTECT(nProtected);
	//Rprintf("Return\n");
	return(ANS);
}

