/* version  5/20/2005 */
/* Rsample.c progrma for sampling without replacement in R  MC 11/2002 */
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

void   update_MCMC_freq(double *MCMC_probs, int *model, int p, int m);
double cond_prob(double *model, int j, int n, double *mean, double *beta_matrix , double eps);

void update_cond_tree(SEXP modelspace, struct Node *tree, SEXP modeldim, struct Var *vars, int p, int n, int kt, int *model, double *real_model, double *marg_probs, double *beta_matrix, double eps);

void  update_Cov(double *Cov, double *priorCov, double *SSgam, double *marg_probs, int n, int m, int print);

void insert_model_tree(struct Node *tree, struct Var *vars,  int n, int *model, int num_models);

// [[register]]
SEXP mcmc(SEXP Y, SEXP X, SEXP Rprobinit, SEXP Rmodeldim, SEXP incint, SEXP Ralpha,
          SEXP method, SEXP modelprior, SEXP Rupdate, SEXP Rbestmodel, SEXP plocal,
          SEXP BURNIN_Iterations, SEXP MCMC_Iterations, SEXP LAMBDA, SEXP DELTA, SEXP Rthin)
{
  SEXP   Rse_m, Rcoef_m, Rmodel_m;

  SEXP   RXwork = PROTECT(duplicate(X)), RYwork = PROTECT(duplicate(Y));
  SEXP   Rbestmarg = PROTECT(allocVector(REALSXP, 1));
  int nProtected = 3, nUnique=0, newmodel=0;
  int nModels=LENGTH(Rmodeldim);

  //  Rprintf("Allocating Space for %d Models\n", nModels) ;
  SEXP ANS = PROTECT(allocVector(VECSXP, 15)); ++nProtected;
  SEXP ANS_names = PROTECT(allocVector(STRSXP, 15)); ++nProtected;
  SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
  SEXP MCMCprobs= PROTECT(duplicate(Rprobinit)); ++nProtected;
  SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
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

  double *Xwork, *Ywork, *coefficients,*probs, shrinkage_m, *MCMC_probs,
    SSY, yty, ybar, mse_m, *se_m, MH=0.0, prior_m=1.0, *real_model,
    R2_m, RSquareFull, alpha, logmargy, postold, postnew;
  int nobs, p, k, i, j, m, n, l, pmodel, pmodel_old, *xdims, *model_m, *bestmodel, *varin, *varout;
  int mcurrent,  update, n_sure;
  double  problocal, *pigamma,  eps, *hyper_parameters;
  double *XtX, *XtY, *XtXwork, *XtYwork, *SSgam, *Cov, *priorCov, *marg_probs;
  double one=1.0, zero=0.0, lambda,  delta;

  int inc=1, p2, thin;
  int *model, *modelold, bit, *modelwork, old_loc, new_loc;
  char uplo[] = "U", trans[]="T";
  struct Var *vars;	/* Info about the model variables. */
  NODEPTR tree, branch;

  /* get dimsensions of all variables */


  nobs = LENGTH(Y);
  xdims = INTEGER(getAttrib(X,R_DimSymbol));
  p = xdims[1];
  k = LENGTH(modelprobs);
  update = INTEGER(Rupdate)[0];
  lambda=REAL(LAMBDA)[0];
  delta = REAL(DELTA)[0];
  thin = INTEGER(Rthin)[0];
  //  Rprintf("delta %f lambda %f", delta, lambda);
  eps = DBL_EPSILON;
  problocal = REAL(plocal)[0];
  //  Rprintf("Update %i and prob.switch %f\n", update, problocal);
  /* Extract prior on models  */
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));

  /*  Rprintf("n %d p %d \n", nobs, p);  */

  Ywork = REAL(RYwork);
  Xwork = REAL(RXwork);


 /* Allocate other variables.  */
  XtX  = (double *) R_alloc(p * p, sizeof(double));
  XtXwork  = (double *) R_alloc(p * p, sizeof(double));
  XtY = vecalloc(p);
  XtYwork = vecalloc(p);


  /* create X matrix */
  for (j=0, l=0; j < p; j++) {
    for (i = 0; i < p; i++) {
      XtX[j*p + i] = 0.0;
    }
    /*    for (i=0; i < nobs; i++) {
       Xmat[i][j] =  REAL(X)[l];
       Xwork[l] = Xmat[i][j];
       l = l + 1;
       } */
  }
  //  PROTECT(Rprobs = NEW_NUMERIC(p));
  //initprobs = REAL(Rprobinit);


 p2 = p*p;
 ybar = 0.0; SSY = 0.0; yty = 0.0;


 F77_NAME(dsyrk)(uplo, trans, &p, &nobs, &one, &Xwork[0], &nobs, &zero, &XtX[0], &p);
 yty = F77_NAME(ddot)(&nobs, &Ywork[0], &inc, &Ywork[0], &inc);
 for (i = 0; i< nobs; i++) {
     ybar += Ywork[i];
  }

  ybar = ybar/ (double) nobs;
  SSY = yty - (double) nobs* ybar *ybar;

  F77_NAME(dgemv)(trans, &nobs, &p, &one, &Xwork[0], &nobs, &Ywork[0], &inc, &zero, &XtY[0],&inc);

  alpha = REAL(Ralpha)[0];

  vars = (struct Var *) R_alloc(p, sizeof(struct Var));
  probs =  REAL(Rprobs);
  n = sortvars(vars, probs, p);

  for (i =n; i <p; i++) REAL(MCMCprobs)[vars[i].index] = probs[vars[i].index];
  for (i =0; i <n; i++) REAL(MCMCprobs)[vars[i].index] = 0.0;
  MCMC_probs =  REAL(MCMCprobs);


  pigamma = vecalloc(p);
  real_model = vecalloc(n);
  marg_probs = vecalloc(n);
  modelold = ivecalloc(p);
  model = ivecalloc(p);
  modelwork= ivecalloc(p);
  varin= ivecalloc(p);
  varout= ivecalloc(p);


  /* create gamma gamma' matrix */
  SSgam  = (double *) R_alloc(n * n, sizeof(double));
  Cov  = (double *) R_alloc(n * n, sizeof(double));
  priorCov  = (double *) R_alloc(n * n, sizeof(double));
  for (j=0; j < n; j++) {
    for (i = 0; i < n; i++) {
      SSgam[j*n + i] = 0.0;
      Cov[j*n + i] = 0.0;
      priorCov[j*n + i] = 0.0;
      if (j == i)  priorCov[j*n + i] = lambda;
    }
    marg_probs[i] = 0.0;
  }


  /* Make space for the models and working variables. */

  /*  pivot = ivecalloc(p);
  qraux = vecalloc(p);
  work =  vecalloc(2 * p);
  effects = vecalloc(nobs);
  v =  vecalloc(p * p);
  betaols = vecalloc(p);
  */



  /*  Rprintf("Fit Full Model\n"); */

  if (nobs <= p) {RSquareFull = 1.0;}
  else {
  PROTECT(Rcoef_m = NEW_NUMERIC(p));
  PROTECT(Rse_m = NEW_NUMERIC(p));
  coefficients = REAL(Rcoef_m);
  se_m = REAL(Rse_m);
  memcpy(coefficients, XtY,  p*sizeof(double));
  memcpy(XtXwork, XtX, p2*sizeof(double));
  memcpy(XtYwork, XtY,  p*sizeof(double));

  mse_m = yty;
  cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, p, nobs);

  /*olsreg(Ywork, Xwork,  coefficients, se_m, &mse_m, &p, &nobs, pivot,qraux,work,residuals,effects,v, betaols); */
  RSquareFull =  1.0 - (mse_m * (double) ( nobs - p))/SSY;
  UNPROTECT(2);
  }


  /* fill in the sure things */
  for (i = n, n_sure = 0; i < p; i++)  {
      model[vars[i].index] = (int) vars[i].prob;
      if (model[vars[i].index] == 1) ++n_sure;
  }


  GetRNGstate();
  tree = make_node(-1.0);

  /*  Rprintf("For m=0, Initialize Tree with initial Model\n");  */

  m = 0;
  bestmodel = INTEGER(Rbestmodel);

  INTEGER(modeldim)[m] = n_sure;

  /* Rprintf("Create Tree\n"); */
   branch = tree;

   for (i = 0; i< n; i++) {
		bit =  bestmodel[vars[i].index];
		if (bit == 1) {
			if (i < n-1 && branch->one == NULL)
				branch->one = make_node(-1.0);
			if (i == n-1 && branch->one == NULL)
				branch->one = make_node(0.0);
			branch = branch->one;
		} else {
			if (i < n-1 && branch->zero == NULL)
				branch->zero = make_node(-1.0);
			if (i == n-1 && branch->zero == NULL)
				branch->zero = make_node(0.0);
			branch = branch->zero;
		}

		model[vars[i].index] = bit;
		INTEGER(modeldim)[m]  += bit;
		branch->where = 0;
	}



    /*    Rprintf("Now get model specific calculations \n"); */

    pmodel = INTEGER(modeldim)[m];
    PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
    model_m = INTEGER(Rmodel_m);

    for (j = 0, l=0; j < p; j++) {
		if (model[j] == 1) {
            model_m[l] = j;
			l +=1;
		}
    }

    SET_ELEMENT(modelspace, m, Rmodel_m);

    Rcoef_m = NEW_NUMERIC(pmodel); PROTECT(Rcoef_m);
    Rse_m = NEW_NUMERIC(pmodel);   PROTECT(Rse_m);
    coefficients = REAL(Rcoef_m);
    se_m = REAL(Rse_m);

    for (j=0, l=0; j < pmodel; j++) {
		XtYwork[j] = XtY[model_m[j]];
        for  ( i = 0; i < pmodel; i++) {
			XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
		}
    }

    mse_m = yty;
    memcpy(coefficients, XtYwork, sizeof(double)*pmodel);
    cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);

    R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;

    SET_ELEMENT(beta, m, Rcoef_m);
    SET_ELEMENT(se, m, Rse_m);

    REAL(R2)[m] = R2_m;
    REAL(mse)[m] = mse_m;

    gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);

    REAL(sampleprobs)[m] = 1.0;
    REAL(logmarg)[m] = logmargy;
    REAL(shrinkage)[m] = shrinkage_m;
    prior_m  = compute_prior_probs(model,pmodel,p, modelprior);
    REAL(priorprobs)[m] = prior_m;
    REAL(Rbestmarg)[0] = REAL(logmarg)[m];
    UNPROTECT(3);


    old_loc = 0;
    pmodel_old = pmodel;
    nUnique=1;
    INTEGER(counts)[0] = 0;
    postold =  REAL(logmarg)[m] + log(REAL(priorprobs)[m]);
    memcpy(modelold, model, sizeof(int)*p);
  /*   Rprintf("model %d max logmarg %lf\n", m, REAL(logmarg)[m]); */

    /*  Rprintf("Now Sample the Rest of the Models \n");  */


  m = 0;

  while (nUnique < k && m < INTEGER(BURNIN_Iterations)[0]) {

    memcpy(model, modelold, sizeof(int)*p);
    pmodel =  n_sure;
    MH = 1.0;

    if (pmodel_old == n_sure || pmodel_old == n_sure + n){
		MH =  random_walk(model, vars,  n);
		MH =  1.0 - problocal;
    } else {
		if (unif_rand() < problocal) {
			// random
			MH =  random_switch(model, vars, n, pmodel_old, varin, varout );
		} else {
			// Randomw walk proposal flip bit//
			MH =  random_walk(model, vars,  n);
		}
    }

    branch = tree;
    newmodel= 0;

    for (i = 0; i< n; i++) {
		bit =  model[vars[i].index];
		if (bit == 1) {
			if (branch->one != NULL) branch = branch->one;
			else newmodel = 1;
		} else {
			if (branch->zero != NULL)  branch = branch->zero;
			else newmodel = 1.0;
		}
		pmodel  += bit;
    }

    if (pmodel  == n_sure || pmodel == n + n_sure)  MH = 1.0/(1.0 - problocal);

    if (newmodel == 1) {
		new_loc = nUnique;
		PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
		model_m = INTEGER(Rmodel_m);
		for (j = 0, l=0; j < p; j++) {
			if (model[j] == 1) {
				model_m[l] = j;
				l +=1;}
		}

		Rcoef_m = NEW_NUMERIC(pmodel); PROTECT(Rcoef_m);
		Rse_m = NEW_NUMERIC(pmodel);   PROTECT(Rse_m);
		coefficients = REAL(Rcoef_m);
		se_m = REAL(Rse_m);
		for (j=0, l=0; j < pmodel; j++) {
			XtYwork[j] = XtY[model_m[j]];
			for  ( i = 0; i < pmodel; i++) {
				XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
			}
		}

		mse_m = yty;
		memcpy(coefficients, XtYwork, sizeof(double)*pmodel);
		cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);

		R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;
		prior_m = compute_prior_probs(model,pmodel,p, modelprior);
		gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
		postnew = logmargy + log(prior_m);
    } else {
		new_loc = branch->where;
		postnew =  REAL(logmarg)[new_loc] + log(REAL(priorprobs)[new_loc]);
    }

    MH *= exp(postnew - postold);
    //    Rprintf("MH new %lf old %lf\n", postnew, postold);
    if (unif_rand() < MH) {
		if (newmodel == 1)  {
			new_loc = nUnique;
			insert_model_tree(tree, vars, n, model, nUnique);

			INTEGER(modeldim)[nUnique] = pmodel;
			SET_ELEMENT(modelspace, nUnique, Rmodel_m);

			SET_ELEMENT(beta, nUnique, Rcoef_m);
			SET_ELEMENT(se, nUnique, Rse_m);

			REAL(R2)[nUnique] = R2_m;
			REAL(mse)[nUnique] = mse_m;
			REAL(sampleprobs)[nUnique] = 1.0;
			REAL(logmarg)[nUnique] = logmargy;
			REAL(shrinkage)[nUnique] = shrinkage_m;
			REAL(priorprobs)[nUnique] = prior_m;
			UNPROTECT(3);
			++nUnique;
		}

		old_loc = new_loc;
		postold = postnew;
		pmodel_old = pmodel;
		memcpy(modelold, model, sizeof(int)*p);
    } else  {
		if (newmodel == 1) UNPROTECT(3);
    }

    INTEGER(counts)[old_loc] += 1;

    for (i = 0; i < n; i++) {
		/* store in opposite order so nth variable is first */
		real_model[n-1-i] = (double) modelold[vars[i].index];
		REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
	}

	// Update SSgam = gamma gamma^T + SSgam
    F77_NAME(dsyr)("U", &n,  &one, &real_model[0], &inc,  &SSgam[0], &n);
    m++;
  }

  for (i = 0; i < n; i++) {
     REAL(MCMCprobs)[vars[i].index] /= (double) m;
  }
  //  Rprintf("\n%d \n", nUnique);


  // Compute marginal probabilities
  mcurrent = nUnique;
  compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
  compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);



  //  Now sample W/O Replacement
  //  Rprintf("NumUnique Models Accepted %d \n", nUnique);
  INTEGER(NumUnique)[0] = nUnique;


	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));

	if (nUnique < nModels) {
		SEXP modelspaceP = PROTECT(allocVector(VECSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			SEXP model_temp = PROTECT(VECTOR_ELT(modelspace, i));
			SET_ELEMENT(modelspaceP, i, model_temp);
			UNPROTECT(1);
		}
		SET_VECTOR_ELT(ANS, 1, modelspaceP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 1, modelspace);
	}
	SET_STRING_ELT(ANS_names, 1, mkChar("which"));

	if (nUnique < nModels) {
		SEXP logmargP = PROTECT(allocVector(REALSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			REAL(logmargP)[i] = REAL(logmarg)[i];
		}
		SET_VECTOR_ELT(ANS, 2, logmargP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 2, logmarg);
	}
	SET_STRING_ELT(ANS_names, 2, mkChar("logmarg"));

	if (nUnique < nModels) {
		SEXP modelprobsP = PROTECT(allocVector(REALSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			REAL(modelprobsP)[i] = REAL(modelprobs)[i];
		}
		SET_VECTOR_ELT(ANS, 3, modelprobsP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 3, modelprobs);
	}
	SET_STRING_ELT(ANS_names, 3, mkChar("postprobs"));

	if (nUnique < nModels) {
		SEXP priorprobsP = PROTECT(allocVector(REALSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			REAL(priorprobsP)[i] = REAL(priorprobs)[i];
		}
		SET_VECTOR_ELT(ANS, 4, priorprobsP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 4, priorprobs);
	}
	SET_STRING_ELT(ANS_names, 4, mkChar("priorprobs"));

	if (nUnique < nModels) {
		SEXP sampleprobsP = PROTECT(allocVector(REALSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			REAL(sampleprobsP)[i] = REAL(sampleprobs)[i];
		}
		SET_VECTOR_ELT(ANS, 5, sampleprobsP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 5, sampleprobs);
	}
	SET_STRING_ELT(ANS_names, 5, mkChar("sampleprobs"));

	if (nUnique < nModels) {
		SEXP mseP = PROTECT(allocVector(REALSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			REAL(mseP)[i] = REAL(mse)[i];
		}
		SET_VECTOR_ELT(ANS, 6, mseP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 6, mse);
	}
	SET_STRING_ELT(ANS_names, 6, mkChar("mse"));

	if (nUnique < nModels) {
		SEXP betaP = PROTECT(allocVector(VECSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			SEXP beta_temp = PROTECT(VECTOR_ELT(beta, i));
			SET_ELEMENT(betaP, i, beta_temp);
			UNPROTECT(1);
		}
		SET_VECTOR_ELT(ANS, 7, betaP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 7, beta);
	}
	SET_STRING_ELT(ANS_names, 7, mkChar("mle"));

	if (nUnique < nModels) {
		SEXP seP = PROTECT(allocVector(VECSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			SEXP se_temp = PROTECT(VECTOR_ELT(se, i));
			SET_ELEMENT(seP, i, se_temp);
			UNPROTECT(1);
		}
		SET_VECTOR_ELT(ANS, 8, seP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 8, se);
	}
	SET_STRING_ELT(ANS_names, 8, mkChar("mle.se"));

	if (nUnique < nModels) {
		SEXP shrinkageP = PROTECT(allocVector(REALSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			REAL(shrinkageP)[i] = REAL(shrinkage)[i];
		}
		SET_VECTOR_ELT(ANS, 9, shrinkageP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 9, shrinkage);
	}
	SET_STRING_ELT(ANS_names, 9, mkChar("shrinkage"));

	if (nUnique < nModels) {
		SEXP modeldimP = PROTECT(allocVector(INTSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			INTEGER(modeldimP)[i] = INTEGER(modeldim)[i];
		}
		SET_VECTOR_ELT(ANS, 10, modeldimP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 10, modeldim);
	}
	SET_STRING_ELT(ANS_names, 10, mkChar("size"));

	if (nUnique < nModels) {
		SEXP R2P = PROTECT(allocVector(REALSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			REAL(R2P)[i] = REAL(R2)[i];
		}
		SET_VECTOR_ELT(ANS, 11, R2P);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 11, R2);
	}
	SET_STRING_ELT(ANS_names, 11, mkChar("R2"));

	if (nUnique < nModels) {
		SEXP countsP = PROTECT(allocVector(INTSXP, nUnique));
		for (i =0; i < nUnique; i++) {
			INTEGER(countsP)[i] = INTEGER(counts)[i];
		}
		SET_VECTOR_ELT(ANS, 12, countsP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 12, counts);
	}
	SET_STRING_ELT(ANS_names, 12, mkChar("freq"));

	SET_VECTOR_ELT(ANS, 13, MCMCprobs);
	SET_STRING_ELT(ANS_names, 13, mkChar("probne0.MCMC"));

	SET_VECTOR_ELT(ANS, 14, NumUnique);
	SET_STRING_ELT(ANS_names, 14, mkChar("n.Unique"));

	setAttrib(ANS, R_NamesSymbol, ANS_names);
	UNPROTECT(nProtected);
	//	Rprintf("Return\n");
	PutRNGstate();

	return(ANS);
}




