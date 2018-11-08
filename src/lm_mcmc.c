#include "bas.h"

void insert_model_tree(struct Node *tree, struct Var *vars,  int n, int *model, int num_models);

int *GetModel_m(SEXP Rmodel_m, int *model, int p) {
	int *model_m = INTEGER(Rmodel_m);
	for (int j = 0, l=0; j < p; j++) {
		if (model[j] == 1) {
			model_m[l++] = j;
		}
	}
	return model_m;
}

void CreateTree(NODEPTR branch, struct Var *vars, int *bestmodel,
                int *model, int n, int m, SEXP modeldim, SEXP Rparents) {
	for (int i = 0; i< n; i++) {
		int bit =  bestmodel[vars[i].index];
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
}
//caller is responsible to call UNPROTECT(2);
double FitModel(SEXP Rcoef_m, SEXP Rse_m, double *XtY, double *XtX, int *model_m,
			  double *XtYwork, double *XtXwork, double yty, double SSY, int pmodel, int p,
			  int nobs, int m, double *pmse_m, int *rank_m, int pivot, double tol) {

	double *coefficients = REAL(Rcoef_m);
	double *se_m = REAL(Rse_m);


	for (int j=0; j < pmodel; j++) { //subsetting matrix
		XtYwork[j] = XtY[model_m[j]];
		for  ( int i = 0; i < pmodel; i++) {
			XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
		}
	}
	*pmse_m = yty;
	 memcpy(coefficients, XtYwork, sizeof(double)*pmodel);
	if (pivot == 1) {
  	*rank_m = cholregpivot(XtYwork, XtXwork, coefficients, se_m, pmse_m, pmodel, nobs, pivot, tol);
	}
	else {
	  cholreg(XtYwork, XtXwork, coefficients, se_m, pmse_m, pmodel, nobs);
	  *rank_m = pmodel;
}

	 double R2_m = 1.0 - (*pmse_m * (double) ( nobs - *rank_m))/SSY;
//	 Rprintf("R2 = %lf, mse = %lf, rank = %d\n", R2_m, *pmse_m, *rank_m);
        if (R2_m < 0.0 || *rank_m == 1) R2_m = 0.0;

	return R2_m;
}

void SetModel2(double logmargy, double shrinkage_m, double prior_m,
			  SEXP sampleprobs, SEXP logmarg, SEXP shrinkage, SEXP priorprobs, int m) {
	REAL(sampleprobs)[m] = 1.0;
	REAL(logmarg)[m] = logmargy;
	REAL(shrinkage)[m] = shrinkage_m;
	REAL(priorprobs)[m] = prior_m;
}

void SetModel(SEXP Rcoef_m, SEXP Rse_m, SEXP Rmodel_m, double mse_m, double R2_m,
              SEXP beta, SEXP se, SEXP modelspace, SEXP mse, SEXP R2, int m) {

  SET_ELEMENT(beta, m, Rcoef_m);
  SET_ELEMENT(se, m, Rse_m);
  SET_ELEMENT(modelspace, m, Rmodel_m);

  REAL(R2)[m] = R2_m;
  REAL(mse)[m] = mse_m;

  UNPROTECT(3);
}

void SetModel_lm(SEXP Rcoef_m, SEXP Rse_m, SEXP Rmodel_m, double mse_m, double R2_m,
              SEXP beta, SEXP se, SEXP modelspace, SEXP mse, SEXP R2, int m) {

  SET_ELEMENT(beta, m, Rcoef_m);
  SET_ELEMENT(se, m, Rse_m);
  SET_ELEMENT(modelspace, m, Rmodel_m);

  REAL(R2)[m] = R2_m;
  REAL(mse)[m] = mse_m;

}


double GetNextModelCandidate(int pmodel_old, int n, int n_sure, int *model, struct Var *vars, double problocal,
							 int *varin, int *varout, SEXP Rparents) {
	double MH = 1.0;
	if (pmodel_old == n_sure || pmodel_old == n_sure + n){
		MH =  random_walk_heredity(model, vars,  n, Rparents);
		MH =  1.0 - problocal;
	} else {
		if (unif_rand() < problocal) {
			// random
			MH =  random_switch_heredity(model, vars, n, pmodel_old, varin, varout, Rparents );
		} else {
			// Random walk proposal flip bit//
			MH =  random_walk_heredity(model, vars,  n, Rparents);
		}
	}
	return MH;
}

double random_walk(int *model, struct Var *vars, int n) {
  int index;
  index = ftrunc(n*unif_rand());
  model[vars[index].index] = 1 - model[vars[index].index];
  return(1.0);
}

double random_switch(int *model, struct Var *vars, int n, int pmodel, int *varin, int *varout) {
  int  j, k, swapin, swapout, num_to_swap_in, num_to_swap_out;


  j = 0; k = 0;
  while (j < n && k < pmodel)
  {
    if (model[vars[j].index]==1) {varin[k] = vars[j].index; k++;}
    j++ ;
  }
  num_to_swap_in = k;
  j = 0; k = 0;

  while (j< n)
  {
    if (model[vars[j].index]==0) {varout[k] = vars[j].index; k++;}
    j++ ;
  }
  num_to_swap_out = k;

  swapin = ftrunc(unif_rand()*num_to_swap_in);    // swapin :corresponds to position of randomly chosen included variable
  swapout = ftrunc(unif_rand()*num_to_swap_out);  // swapout :corresponds to position of randomly chosen excluded variable

  model[varin[swapin]] = 0;
  model[varout[swapout]] =1;


  return(1.0);
}

double random_walk_heredity(int *model, struct Var *vars, int n, SEXP Rparents) {
  int index,p,j;
  double *parents;

  parents = REAL(Rparents);

  index = ftrunc(n*unif_rand());
  model[vars[index].index] = 1 - model[vars[index].index];


 int *dims = INTEGER(getAttrib(Rparents,R_DimSymbol));
 p = dims[0];


  if (p > 1) {
   // force in parents
 //  Rprintf("%d %d %d %d\n",n,p,  vars[index].index, model[vars[index].index]);
    if (model[vars[index].index] == 1) {
//  traverse row index of parents to add any missing parents/sibs
    for (j = 0; j < p; j++) {
 //     Rprintf("%d ", (int) parents[vars[index].index*p + j]);
      if (parents[vars[index].index + p*j] == 1.0) {
         model[j] =  model[vars[index].index];
      }
    }}
    else {
//  to drop index, traverse column of parents to identify children/sibs
//  that also need to be dropped.
      for (j = 0; j < p; j++) {
  //      Rprintf("%d ", (int) parents[vars[index].index + p*j]);
        if (parents[vars[index].index*p +j] == 1.0) {
          model[j] =  model[vars[index].index];
        } }
    }
//    Rprintf("\n");

  }
  return(1.0);
}

double random_switch_heredity(int *model, struct Var *vars, int n,
                              int pmodel, int *varin, int *varout, SEXP Rparents)
  {
  int  j, k, p, swapin, swapout, num_to_swap_in, num_to_swap_out;
  double *parents;

  j = 0; k = 0;
  while (j < n && k < pmodel)
  {
    if (model[vars[j].index]==1) {
      varin[k] = vars[j].index;
      k++;}
    j++ ;
  }
  num_to_swap_in = k;
  j = 0; k = 0;

  while (j< n)
  {
    if (model[vars[j].index]==0) {
      varout[k] = vars[j].index;
      k++;}
    j++ ;
  }
  num_to_swap_out = k;

  swapin = ftrunc(unif_rand()*num_to_swap_in);    // swapin :corresponds to position of randomly chosen included variable
  swapout = ftrunc(unif_rand()*num_to_swap_out);  // swapout :corresponds to position of randomly chosen excluded variable

  model[varin[swapin]] = 0;
  model[varout[swapout]] = 1;

  parents = REAL(Rparents);
  int *dims = INTEGER(getAttrib(Rparents,R_DimSymbol));
  p = dims[0];

  // force in parents and sibs of variable that was swapped in

  if (p > 1) {
    //  to drop swapin, traverse column of parents to identify children/sibs
    //  that also need to be dropped.  ignore others
      for (j = 0; j < p; j++) {
        if (parents[varin[swapin]*p +j] == 1.0)   model[j] = 0;
        }

    //  now traverse row of added variable in parents to add any missing parents/sibs
      for (j = 0; j < p; j++) {
        if (parents[varout[swapout] + p*j] == 1.0)   model[j] = 1;
      }
  }
  return(1.0);
  }


// [[register]]
extern SEXP mcmc_new(SEXP Y, SEXP X, SEXP Rweights, SEXP Rprobinit, SEXP Rmodeldim, SEXP incint,
                     SEXP Ralpha,SEXP method,SEXP modelprior, SEXP Rupdate,
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

	//  Rprintf("Allocating Space for %d Models\n", nModels) ;
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

	double *Xwork, *Ywork,*wts, *probs, shrinkage_m,
		mse_m, MH=0.0, prior_m=1.0,
		R2_m, RSquareFull, logmargy, postold=0.0, postnew=0.0;
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



	//	Rprintf("delta %f lambda %f", delta, lambda);

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


	// fill in the sure things
	int *model = ivecalloc(p);
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
	PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
	PROTECT(Rse_m = NEW_NUMERIC(pmodel));

	model_m = GetModel_m(Rmodel_m, model, p);
	//evaluate logmargy and shrinkage

	R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel,
                 p, nobs, m, &mse_m, &rank_m, pivot, tol);
	INTEGER(rank)[0] = rank_m;

	gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy,
               &shrinkage_m);


	prior_m  = compute_prior_probs(model,pmodel,p, modelprior);
	if (prior_m == 0.0)  Rprintf("warning initial model has 0 prior probabilty\n");
	SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
	SetModel(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2, m);

	int nUnique=0, newmodel=0, nsamples=0;
	double *real_model = vecalloc(n);
	int *modelold = ivecalloc(p);
	int old_loc = 0, new_loc = 0;
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
		  prior_m = compute_prior_probs(model,pmodel,p, modelprior);
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

