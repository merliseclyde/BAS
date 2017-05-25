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
#include "sampling.h"
void update_tree(SEXP modelspace, struct Node *tree, SEXP modeldim, struct Var *vars, int k, int p, int n, int kt, int *model);

double CalculateRSquareFull(double *XtY, double *XtX, double *XtXwork, double *XtYwork,
							SEXP Rcoef_m, SEXP Rse_m, int p, int nobs, double yty, double SSY);
int *GetModel_m(SEXP Rmodel_m, int *model, int p);
double FitModel(SEXP Rcoef_m, SEXP Rse_m, double *XtY, double *XtX, int *model_m,
			  double *XtYwork, double *XtXwork, double yty, double SSY, int pmodel, int p,
			  int nobs, int m, double *pmse_m);
void SetModel2(double logmargy, double shrinkage_m, double prior_m,
			  SEXP sampleprobs, SEXP logmarg, SEXP shrinkage, SEXP priorprobs, int m);
void SetModel(SEXP Rcoef_m, SEXP Rse_m, SEXP Rmodel_m, double mse_m, double R2_m,
			  SEXP beta, SEXP se, SEXP modelspace, SEXP mse, SEXP R2, int m);

void CreateTree_with_pigamma(NODEPTR branch, struct Var *vars, int *bestmodel, int *model, int n, int m, SEXP modeldim, double *pigamma) {
	for (int i = 0; i< n; i++) {
		pigamma[i] = 1.0;
		int bit =  bestmodel[vars[i].index];
		if (bit == 1) {
			for (int j=0; j<=i; j++)  pigamma[j] *= branch->prob;
			if (i < n-1 && branch->one == NULL)
				branch->one = make_node(vars[i+1].prob);
			if (i == n-1 && branch->one == NULL)
				branch->one = make_node(0.0);
			branch = branch->one;
		} else {
			for (int j=0; j<=i; j++)  pigamma[j] *= (1.0 - branch->prob);
			if (i < n-1 && branch->zero == NULL)
				branch->zero = make_node(vars[i+1].prob);
			if (i == n-1 && branch->zero == NULL)
				branch->zero = make_node(0.0);
			branch = branch->zero;
		}
		model[vars[i].index] = bit;
		INTEGER(modeldim)[m]  += bit;
	}
}
void Substract_visited_probability_mass(NODEPTR branch, struct Var *vars, int *model, int n, int m,  double *pigamma, double eps) {
	for (int i = 0; i < n; i++) {
		int bit = model[vars[i].index];
		double prone = branch->prob;
		if (bit == 1) prone -= pigamma[i];
		double denom = 1.0 - pigamma[i];
		if (denom <= 0.0) {
			if (denom < 0.0) {
				Rprintf("neg denominator %le %le %le !!!\n", pigamma, denom, prone);
				if (branch->prob < 0.0 && branch->prob < 1.0) {
					Rprintf("non extreme %le\n", branch->prob);
				}
			}
			denom = 0.0;
		} else {
			if  (prone <= 0)  prone = 0.0;
			if  (prone > denom)  {
				if (prone <= eps) prone = 0.0;
				else prone = 1.0;
				/* Rprintf("prone > 1 %le %le %le %le !!!\n", pigamma, denom, prone, eps);*/
			}
			else prone = prone/denom;
		}
		if (prone > 1.0 || prone < 0.0)
			Rprintf("%d %d Probability > 1!!! %le %le  %le %le \n",
			m, i, prone, branch->prob, denom, pigamma);
		branch->prob  = prone;
		if (bit == 1) branch = branch->one;
		else  branch = branch->zero;
	}
}
void GetNextModel_swop(NODEPTR branch, struct Var *vars,
                       int *model, int n, int m,  double *pigamma,
                   		 double problocal, SEXP modeldim, int *bestmodel) {
	for (int i = 0; i< n; i++) {
		pigamma[i] = 1.0;
		int bit =  withprob(branch->prob);
		int local = withprob(problocal);
		if ( local == 1 && (branch->prob < .999) && (branch->prob > 0.001)) {
			bit = bestmodel[vars[i].index];
		} else {
			bit =  withprob(branch->prob);
		}
		if (bit == 1) {
			for (int j=0; j<=i; j++)  pigamma[j] *= branch->prob;
			if (i < n-1 && branch->one == NULL) branch->one = make_node(vars[i+1].prob);
			if (i == n-1 && branch->one == NULL) branch->one = make_node(0.0);
			branch = branch->one;
		} else {
			for (int j=0; j<=i; j++)  pigamma[j] *= (1.0 - branch->prob);
			if (i < n-1 && branch->zero == NULL) branch->zero = make_node(vars[i+1].prob);
			if (i == n-1 && branch->zero == NULL) branch->zero = make_node(0.0);
			branch = branch->zero;
		}
		model[vars[i].index] = bit;
		INTEGER(modeldim)[m]  += bit;
	}
}


extern SEXP sampleworep_new(SEXP Y, SEXP X, SEXP Rweights, SEXP Rprobinit, SEXP Rmodeldim, SEXP incint, SEXP Ralpha,SEXP method, SEXP modelprior, SEXP Rupdate, SEXP Rbestmodel, SEXP plocal) {
	int nProtected = 0;
	SEXP RXwork = PROTECT(duplicate(X)); nProtected++;
	SEXP RYwork = PROTECT(duplicate(Y)); nProtected++;
	int nModels=LENGTH(Rmodeldim);

	//  Rprintf("Allocating Space for %d Models\n", nModels) ;
	SEXP ANS = PROTECT(allocVector(VECSXP, 12)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 12)); ++nProtected;
	SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
	SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP modeldim =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
	SEXP beta = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP se = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP mse = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP priorprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP logmarg = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP sampleprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;

	double *Xwork, *Ywork, *wts, *probs, shrinkage_m, mse_m, R2_m, RSquareFull, Rbestmarg, logmargy;
	int i;

	//get dimsensions of all variables
	int nobs = LENGTH(Y);
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int k = LENGTH(modelprobs);
	double alpha = REAL(Ralpha)[0];
	SEXP Rbestmodel_new = PROTECT(duplicate(Rbestmodel)); nProtected++;

	int update = INTEGER(Rupdate)[0];
	double eps = DBL_EPSILON;
	double problocal = REAL(plocal)[0];

//	memset(INTEGER(modeldim), 1, k*sizeof(int));
	Ywork = REAL(RYwork);
	Xwork = REAL(RXwork);
	wts = REAL(Rweights);

	double *XtXwork, *XtYwork,*XtX, *XtY, yty,SSY;
	PrecomputeData(Xwork, Ywork, wts, &XtXwork, &XtYwork, &XtX, &XtY, &yty, &SSY, p, nobs);

	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables.
	probs =  REAL(Rprobs);
	int n = sortvars(vars, probs, p);

	SEXP  Rse_m = NULL, Rcoef_m = NULL;
	RSquareFull = CalculateRSquareFull(XtY, XtX, XtXwork, XtYwork, Rcoef_m, Rse_m, p, nobs, yty, SSY);

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
	int *bestmodel = INTEGER(Rbestmodel_new);
	REAL(logmarg)[m] = 0.0;

	for (i = n; i < p; i++)  {
		model[vars[i].index] = bestmodel[vars[i].index];
		INTEGER(modeldim)[m]  +=  bestmodel[vars[i].index];
	}

	double *pigamma = vecalloc(p);
	branch = tree;
	CreateTree_with_pigamma(branch, vars, bestmodel, model, n, m, modeldim,pigamma);

	branch=tree;
	Substract_visited_probability_mass(branch, vars, model, n, m, pigamma,eps);

	int pmodel = INTEGER(modeldim)[m];

	SEXP Rmodel_m = PROTECT(NEW_INTEGER(pmodel));
	memset(INTEGER(Rmodel_m), 0, pmodel * sizeof(int));
	PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
	PROTECT(Rse_m = NEW_NUMERIC(pmodel));
  GetModel_m(Rmodel_m, model, p);
  int *model_m = INTEGER(Rmodel_m);
	R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel, p, nobs, m, &mse_m);
	gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
	double prior_m  = compute_prior_probs(model,pmodel,p, modelprior);



	SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
	SetModel(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2, m);
	//Rprintf("model %d max logmarg %lf\n", m, REAL(logmarg)[m]);

	Rbestmarg = REAL(logmarg)[m];

	int *modelwork= ivecalloc(p);

	// Sample models
	for (m = 1;  m < k; m++) {
	  INTEGER(modeldim)[m] = 0;
		for (i = n; i < p; i++)  {
			INTEGER(modeldim)[m]  +=  model[vars[i].index];
		}

		branch = tree;
		GetNextModel_swop(branch, vars, model, n, m, pigamma, problocal, modeldim, bestmodel);

		/* Now subtract off the visited probability mass. */
		branch=tree;
		Substract_visited_probability_mass(branch, vars, model, n, m, pigamma,eps);

		/* Now get model specific calculations */
		pmodel = INTEGER(modeldim)[m];
		if (pmodel < 1) Rprintf("%f\n", pmodel);
		SEXP Rmodel_m = PROTECT(NEW_INTEGER(pmodel));
		memset(INTEGER(Rmodel_m), 0, pmodel * sizeof(int));
		PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
		PROTECT(Rse_m = NEW_NUMERIC(pmodel));
		GetModel_m(Rmodel_m, model, p);
    int *model_m = INTEGER(Rmodel_m);

		R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel, p, nobs, m, &mse_m);
		gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
		prior_m = compute_prior_probs(model,pmodel,p, modelprior);
		SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
		SetModel(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2,m);
		REAL(sampleprobs)[m] = pigamma[0];

		//update best model
		if (REAL(logmarg)[m] > Rbestmarg) {
			for (i=0; i < p; i++) {
				bestmodel[i] = model[i];
			}
			Rbestmarg = REAL(logmarg)[m];
		}

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

	setAttrib(ANS, R_NamesSymbol, ANS_names);
	PutRNGstate();

	UNPROTECT(nProtected);


	return(ANS);
}
