#include "bas.h"


// [[register]]
SEXP mcmc_new(SEXP Y, SEXP X, SEXP Rweights, SEXP Rprobinit, SEXP Rmodeldim,
              SEXP incint, SEXP Ralpha, SEXP method, SEXP modelprior, SEXP Rupdate,
              SEXP Rbestmodel, SEXP plocal, SEXP BURNIN_Iterations,
              SEXP MCMC_Iterations, SEXP LAMBDA, SEXP DELTA,
              SEXP Rthin, SEXP Rparents, SEXP Rpivot, SEXP Rtol)
{
	int nProtected = 0;
	SEXP RXwork = PROTECT(duplicate(X)); nProtected++;
	SEXP RYwork = PROTECT(duplicate(Y)); nProtected++;
	int nModels=LENGTH(Rmodeldim);
	int pivot = LOGICAL(Rpivot)[0];
	double tol = REAL(Rtol)[0];

  // Rprintf("Allocating Space for %d Models\n", nModels) ;
	SEXP ANS = PROTECT(allocVector(VECSXP, 16)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 16)); ++nProtected;
	SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
	SEXP MCMCprobs= PROTECT(duplicate(Rprobinit)); ++nProtected;
	SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP rank = PROTECT(allocVector(INTSXP, nModels)); ++nProtected;
	SEXP modeldim =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
	SEXP counts =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
	SEXP beta = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP se = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP mse = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP priorprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP logmarg = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP sampleprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
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
	PROTECT_INDEX rank_idx;
	PROTECT_WITH_INDEX(rank, &rank_idx);
	PROTECT_INDEX beta_idx;
	PROTECT_WITH_INDEX(beta, &beta_idx);
	PROTECT_INDEX se_idx;
	PROTECT_WITH_INDEX(se, &se_idx);
	PROTECT_INDEX mse_idx;
	PROTECT_WITH_INDEX(mse, &mse_idx);
	PROTECT_INDEX modelprobs_idx;
	PROTECT_WITH_INDEX(modelprobs, &modelprobs_idx);
	PROTECT_INDEX priorprobs_idx;
	PROTECT_WITH_INDEX(priorprobs, &priorprobs_idx);
	PROTECT_INDEX logmarg_idx;
	PROTECT_WITH_INDEX(logmarg, &logmarg_idx);
	PROTECT_INDEX sampleprobs_idx;
	PROTECT_WITH_INDEX(sampleprobs, &sampleprobs_idx);
	

	double *Xwork, *Ywork,*wts, *probs, shrinkage_m,
		mse_m, MH=0.0, prior_m=1.0,
		R2_m, RSquareFull, logmargy, postold, postnew;
	int i, m, n, pmodel_old, *model_m, *bestmodel, rank_m;
	int mcurrent, n_sure;


	//get dimsensions of all variables
	int nobs = LENGTH(Y);
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int k = LENGTH(modelprobs);
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
	
	for (i =n; i <p; i++) REAL(MCMCprobs)[vars[i].index] = probs[vars[i].index];
	for (i =0; i <n; i++) REAL(MCMCprobs)[vars[i].index] = 0.0;
	
	int noInclusionIs1 = no_prior_inclusion_is_1(p, probs);


	//  allocate working model and fill in the sure things
	int *model = ivecalloc(p);
	memset(model, 0, p * sizeof(int));

	for (i = n, n_sure = 0; i < p; i++)  {
		model[vars[i].index] = (int) vars[i].prob;
		if (model[vars[i].index] == 1) ++n_sure;
	}

	SEXP Rse_m = NULL, Rcoef_m = NULL, Rmodel_m=NULL;
	RSquareFull = CalculateRSquareFull(XtY, XtX, XtXwork, XtYwork, Rcoef_m, Rse_m, p, nobs, yty, SSY);

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
	//evaluate logmargy and shrinkage

	R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel,
                 p, nobs, m, &mse_m, &rank_m, pivot, tol);
	INTEGER(rank)[0] = rank_m;

	gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy,
               &shrinkage_m);


	prior_m  = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
	if (prior_m == 0.0)  error("initial model has 0 prior probabilty\n");
	SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
	SetModel(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2, m);

	int nUnique=0, newmodel=0, nsamples=0;
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
	
	
	while (nUnique < k && m < (INTEGER(MCMC_Iterations)[0] + INTEGER(BURNIN_Iterations)[0])) {

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
		    gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);

		    postnew = logmargy + log(prior_m);
		    MH *= exp(postnew - postold);
		  }}
		else {
		  new_loc = branch->where;
		  postnew =  REAL(logmarg)[new_loc] +
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

		    //record model data
		    SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, nUnique);
		    SetModel(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2,nUnique);

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

		  INTEGER(counts)[old_loc] += 1;

		  for (i = 0; i < n; i++) {
		    // store in opposite order so nth variable is first
		    real_model[n-1-i] = (double) modelold[vars[i].index];
		    REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
		  }
		  nsamples++;
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
	compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);

	INTEGER(NumUnique)[0] = nUnique;

	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));

	if (nUnique < nModels) {
	  nModels = nUnique;
	  SETLENGTH(counts, nUnique);
	  REPROTECT(counts= Rf_lengthgets(counts, nUnique), counts_idx);
	  REPROTECT(logmarg= Rf_lengthgets(logmarg, nUnique), logmarg_idx);
	  REPROTECT(modelprobs= Rf_lengthgets(modelprobs, nUnique), modelprobs_idx);
	  REPROTECT(priorprobs= Rf_lengthgets(priorprobs, nUnique), priorprobs_idx);
	  REPROTECT(sampleprobs= Rf_lengthgets(sampleprobs, nUnique), sampleprobs_idx);
	  REPROTECT(mse = Rf_lengthgets(mse, nUnique), mse_idx);
	  REPROTECT(shrinkage = Rf_lengthgets(shrinkage, nUnique), shrinkage_idx);
	  REPROTECT(modeldim= Rf_lengthgets(modeldim, nUnique), modeldim_idx);
	  REPROTECT(R2= Rf_lengthgets(R2, nUnique), R2_idx);
	  REPROTECT(se= Rf_lengthgets(se, nUnique), se_idx);
	  REPROTECT(rank = Rf_lengthgets(rank, nUnique), rank_idx);
	  REPROTECT(modelspace = Rf_lengthgets(modelspace, nUnique), modelspace_idx);
	  REPROTECT(beta = Rf_lengthgets(beta, nUnique), beta_idx);
	  REPROTECT(se= Rf_lengthgets(se, nUnique), se_idx);
	  
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

SET_VECTOR_ELT(ANS, 13, counts);
SET_STRING_ELT(ANS_names, 13, mkChar("freq"));

SET_VECTOR_ELT(ANS, 14, MCMCprobs);
SET_STRING_ELT(ANS_names, 14, mkChar("probne0.MCMC"));

SET_VECTOR_ELT(ANS, 15, NumUnique);
SET_STRING_ELT(ANS_names, 15, mkChar("n.Unique"));

setAttrib(ANS, R_NamesSymbol, ANS_names);

	PutRNGstate();
    UNPROTECT(nProtected);
    //	Rprintf("Return\n");
	return(ANS);
}

