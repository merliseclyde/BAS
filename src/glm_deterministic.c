#include "sampling.h"
int *GetModel_m(SEXP Rmodel_m, int *model, int p);
double FitModel(SEXP Rcoef_m, SEXP Rse_m, double *XtY, double *XtX, int *model_m,
			  double *XtYwork, double *XtXwork, double yty, double SSY, int pmodel, int p,
			  int nobs, int m, double *pmse_m);
SEXP gglm_lpy(SEXP RX, SEXP RY,SEXP Ra, SEXP Rb, SEXP Rs, SEXP Rcoef, SEXP Rmu);

void SetModel2(double logmargy, double shrinkage_m, double prior_m,
			  SEXP sampleprobs, SEXP logmarg, SEXP shrinkage, SEXP priorprobs, int m);
void SetModel1(SEXP Rfit, SEXP Rmodel_m, 
			  SEXP beta, SEXP se, SEXP modelspace, SEXP deviance, SEXP R2, SEXP Q, int m);
SEXP glm_FitModel(SEXP RX, SEXP RY, SEXP Rmodel_m,  //input data
			  SEXP Roffset, SEXP Rweights, SEXP family, SEXP Rcontrol,
			  SEXP Ra, SEXP Rb, SEXP Rs);

int topk(Bit **models, double *prob, int k, struct Var *vars, int n, int p);
void insert_children(int subset, double *list, double *subsetsum,
					 int *queue, int *queuesize, int *tablesize,
					 int *parent, int *pattern, int *position,
					 int *type, char *bits, int  n);
void do_insert(int child, double *subsetsum, int *queue);
int get_next(double *subsetsum, int *queue, int *queuesize);
void set_bits(char *bits, int subset, int *pattern, int *position, int n);
void print_subset(int subset, int rank, Bit **models, Bit *model,
				  double *subsetsum, int *pattern, int *position, 
				  int n, struct Var *vars, int p);
int withprob(double p);

SEXP glm_deterministic(SEXP Y, SEXP X, SEXP Roffset, SEXP Rweights, 
			  SEXP Rprobinit, SEXP Rmodeldim, SEXP modelprior, 
			  SEXP Ra, SEXP Rb, SEXP Rs,
			  SEXP family, SEXP Rcontrol) {
	int nProtected = 0;
	int nModels=LENGTH(Rmodeldim);

	//  Rprintf("Allocating Space for %d Models\n", nModels) ;
	SEXP ANS = PROTECT(allocVector(VECSXP, 13)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 13)); ++nProtected;
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

	double *probs,logmargy;

	//get dimsensions of all variables 
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int k = LENGTH(modelprobs);

	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables. 
	probs =  REAL(Rprobs);
	int n = sortvars(vars, probs, p); 

	Bit **models = cmatalloc(k,p);
	int *model = (int *) R_alloc(p, sizeof(int));
	k = topk(models, probs, k, vars, n, p);

	/* now fit all top k models */
	for (int m=0; m < k; m++) {
		int pmodel = 0; 
		double pigamma = 1.0;
		for (int j = 0; j < p; j++) { 
			model[j] = (int) models[m][j];
			pmodel += (int) models[m][j];
			pigamma *= (double)((int) models[m][j])*probs[j] + 
				(1.0 - (double)((int) models[m][j]))*(1.0 -  probs[j]);
		}

		SEXP Rmodel_m =	PROTECT(allocVector(INTSXP,pmodel));
		GetModel_m(Rmodel_m, model, p);
		//evaluate logmargy and shrinkage
		SEXP glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights, family, Rcontrol, Ra, Rb, Rs));	
		double prior_m  = compute_prior_probs(model,pmodel,p, modelprior);
		logmargy = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
		SetModel2(logmargy, NA_REAL, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
		REAL(sampleprobs)[m] = pigamma;
		SetModel1(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q,m);
		UNPROTECT(2);
	}

	compute_modelprobs(modelprobs, logmarg, priorprobs, k);
	compute_margprobs_old(models, modelprobs, probs, k, p); 

	/*    freechmat(models,k); */
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
	SET_STRING_ELT(ANS_names, 7, mkChar("coefficients"));

	SET_VECTOR_ELT(ANS, 8, se);
	SET_STRING_ELT(ANS_names, 8, mkChar("se"));

	SET_VECTOR_ELT(ANS, 9, shrinkage);
	SET_STRING_ELT(ANS_names, 9, mkChar("shrinkage"));

	SET_VECTOR_ELT(ANS, 10, modeldim);
	SET_STRING_ELT(ANS_names, 10, mkChar("size"));

	SET_VECTOR_ELT(ANS, 11, R2);
	SET_STRING_ELT(ANS_names, 11, mkChar("R2"));

	SET_VECTOR_ELT(ANS, 12, Q);
	SET_STRING_ELT(ANS_names, 12, mkChar("Q"));

	setAttrib(ANS, R_NamesSymbol, ANS_names);
	UNPROTECT(nProtected);

	return(ANS);  

}
