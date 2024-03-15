#include "bas.h"


// [[register]]
SEXP amcmc(SEXP Y, SEXP X, SEXP Rweights, SEXP Rprobinit, SEXP Rmodeldim,
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

	Rprintf("Allocating Space for %d Models AMCMC\n", nModels) ;
	double *Xwork, *Ywork,*wts, *probs, shrinkage_m,
		mse_m, MH=0.0, prior_m=1.0,
		R2_m, RSquareFull, logmargy, postold, postnew;

	double *SSgam, *Cov, *priorCov, *marg_probs, one=1.0, wt = 1.0;
	int i,j, m, n, pmodel_old, *model_m, *bestmodel, rank_m;
	int mcurrent, n_sure, inc=1;
	int print = 1;

	Rprintf("Allocating Space for %d Models 2 \n", nModels) ;
	Rprintf("AMCMC\n") ;
	
	//get dimsensions of all variables
	int nobs = LENGTH(Y);
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int k = LENGTH(modelprobs);
	double lambda=REAL(LAMBDA)[0];
	double delta = REAL(DELTA)[0];
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

	
	/* create gamma gamma' matrix */
	SSgam  = (double *) R_alloc(n * n, sizeof(double));
	Cov  = (double *) R_alloc(n * n, sizeof(double));
	priorCov  = (double *) R_alloc(n * n, sizeof(double));
	marg_probs = (double *) R_alloc(n, sizeof(double));
	
	memset(SSgam, 0.0, n*n*sizeof(double));
	memset(Cov, 0.0, n*n*sizeof(double));
	memset(priorCov, 0.0, n*n*sizeof(double));
	memset(marg_probs, 0.0, n*sizeof(double));
	
	for (j=0; j < n; j++) {
	  for (i = 0; i < n; i++) {
	    if (j == i)  priorCov[j*n + i] = lambda;
	  }
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
	
	Rprintf("using AMCMC sampling - initialize\n");
	
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
		    if ((m % thin) == 0 )  {

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
//		    real_model[n-1-i] = (double) modelold[vars[i].index];
        real_model[i] = (double) modelold[vars[i].index];
		    REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
		  }
		  
		  // Update SSgam = gamma gamma^T + SSgam 
		  // use Upper Triangular
		  F77_NAME(dsyr)("U", &n,  &one, &real_model[0], &inc,  &SSgam[0], &n FCONE);
		  nsamples++;
		}
		m++;
	}
	
	// Compute marginal probabilities  
	mcurrent = nUnique;
	compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);        
	for (i = 0; i < n; i++) {
	  marg_probs[i] = wt*(REAL(MCMCprobs)[vars[i].index]/ (double) m) + 
	    (1.0 - wt)* probs[vars[i].index];
	}	
	print=0;
	update_Cov(Cov, priorCov, SSgam, marg_probs, n, m, print);
	
	// Global-Proposal
	// Initialize post old proposal
	double pigammaold = 0.0;
	for (i = 0; i < n; i++) {
	  if (modelold[vars[i].index] == 1 ){	
	    real_model[i] = 1.0;
	    pigammaold += log(cond_prob(real_model,i, n, marg_probs,Cov, delta));
	  }
	  else {
	    real_model[i] = 0.0;
	    pigammaold += log(1.0 - cond_prob(real_model,i, n, marg_probs,Cov, delta));
	  }
	}
	
  // now use AMCMC
  
  Rprintf("Now start AMCMC\n");
  
  
  while (nUnique < k && m < INTEGER(MCMC_Iterations)[0]) {
    
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
        if ((m % thin) == 0 )  {
          
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
	  SETLENGTH(modelspace, nUnique);
	  SETLENGTH(logmarg, nUnique);
	  SETLENGTH(modelprobs, nUnique);
	  SETLENGTH(priorprobs, nUnique);
	  SETLENGTH(sampleprobs, nUnique);
	  SETLENGTH(counts, nUnique);
	  SETLENGTH(beta, nUnique);
	  SETLENGTH(se, nUnique);
	  SETLENGTH(mse, nUnique);
	  SETLENGTH(shrinkage, nUnique);
	  SETLENGTH(modeldim, nUnique);
	  SETLENGTH(R2, nUnique);
	  SETLENGTH(rank, nUnique);
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

