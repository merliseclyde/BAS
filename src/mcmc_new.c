#include "sampling.h"

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
double CalculateRSquareFull(double *XtY, double *XtX, double *XtXwork, double *XtYwork, 
							SEXP Rcoef_m, SEXP Rse_m, int p, int nobs, double yty, double SSY) {
	double RSquareFull;
	if (nobs <= p) {
		RSquareFull = 1.0;
	} else {
		PROTECT(Rcoef_m = NEW_NUMERIC(p));
		PROTECT(Rse_m = NEW_NUMERIC(p));
		double *coefficients = REAL(Rcoef_m);  
		double *se_m = REAL(Rse_m);
		memcpy(coefficients, XtY,  p*sizeof(double));
		memcpy(XtXwork, XtX, p * p *sizeof(double));
		memcpy(XtYwork, XtY,  p*sizeof(double));

		double mse_m = yty; 
		cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, p, nobs);  

		RSquareFull =  1.0 - (mse_m * (double) ( nobs - p))/SSY;
		UNPROTECT(2);
	}
	return RSquareFull;
}
void CreateTree(NODEPTR branch, struct Var *vars, int *bestmodel, int *model, int n, int m, SEXP modeldim) {
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
			  int nobs, int m, double *pmse_m) {

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
	cholreg(XtYwork, XtXwork, coefficients, se_m, pmse_m, pmodel, nobs);  

	double R2_m = 1.0 - (*pmse_m * (double) ( nobs - pmodel))/SSY;

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

void PrecomputeData(double *Xwork, double *Ywork, double **pXtXwork, double **pXtYwork, double **pXtX, double **pXtY, double *yty, double *SSY, int p, int nobs) {
	char uplo[] = "U", trans[]="T";
	double one=1.0, zero=0.0;
	int inc=1;
	
	int p2 = p * p;
	*pXtXwork  = (double *) R_alloc(p2, sizeof(double));
	*pXtYwork = vecalloc(p);
	*pXtX  = (double *) R_alloc(p2, sizeof(double));
	*pXtY = vecalloc(p);

	//precompute XtX
	memset(*pXtX,0, p2 * sizeof(double));
	//F77_NAME(dsyrk)(uplo, trans, &p, &nobs, &one, &Xwork[0], &nobs, &zero, &(*pXtX)[0], &p); 
	F77_NAME(dsyrk)(uplo, trans, &p, &nobs, &one, &Xwork[0], &nobs, &zero, *pXtX, &p); 
	*yty = F77_NAME(ddot)(&nobs, &Ywork[0], &inc, &Ywork[0], &inc);
	double ybar = 0.0; 
	for (int i = 0; i< nobs; i++) {
		ybar += Ywork[i];
	}
	ybar = ybar/ (double) nobs;
	*SSY = *yty - (double) nobs* ybar *ybar;
	F77_NAME(dgemv)(trans, &nobs, &p, &one, &Xwork[0], &nobs, &Ywork[0], &inc, &zero, *pXtY,&inc);
}

double GetNextModelCandidate(int pmodel_old, int n, int n_sure, int *model, struct Var *vars, double problocal,
							 int *varin, int *varout) {
	double MH = 1.0;
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
	return MH;
}
SEXP mcmc_new(SEXP Y, SEXP X, SEXP Rprobinit, SEXP Rmodeldim, SEXP incint, SEXP Ralpha,SEXP method, 
			  SEXP modelprior, SEXP Rupdate, SEXP Rbestmodel, SEXP Rbestmarg, SEXP plocal, 
			  SEXP BURNIN_Iterations, SEXP MCMC_Iterations, SEXP LAMBDA, SEXP DELTA)
{
	int nProtected = 0;
	SEXP RXwork = PROTECT(duplicate(X)); nProtected++;
	SEXP RYwork = PROTECT(duplicate(Y)); nProtected++;

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

	double *Xwork, *Ywork, *probs, shrinkage_m,
		mse_m, MH=0.0, prior_m=1.0, 
		R2_m, RSquareFull, logmargy, postold, postnew;
	int i, j, m, n, l, pmodel_old, *model_m, *bestmodel;
	int mcurrent, n_sure;
	
	//get dimsensions of all variables 
	int nobs = LENGTH(Y);
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int k = LENGTH(modelprobs);
	double lambda=REAL(LAMBDA)[0];
	double delta = REAL(DELTA)[0];
	double alpha = REAL(Ralpha)[0];

	//	Rprintf("delta %f lambda %f", delta, lambda);

	Ywork = REAL(RYwork);
	Xwork = REAL(RXwork);

	double *XtXwork, *XtYwork,*XtX, *XtY, yty,SSY;
	PrecomputeData(Xwork, Ywork, &XtXwork, &XtYwork, &XtX, &XtY, &yty, &SSY, p, nobs);

	
	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables. 
	probs =  REAL(Rprobs);
	n = sortvars(vars, probs, p); 
	for (i =n; i <p; i++) REAL(MCMCprobs)[vars[i].index] = probs[vars[i].index];
	for (i =0; i <n; i++) REAL(MCMCprobs)[vars[i].index] = 0.0;
	
	SEXP Rse_m = NULL, Rcoef_m = NULL; 
	RSquareFull = CalculateRSquareFull(XtY, XtX, XtXwork, XtYwork, Rcoef_m, Rse_m, p, nobs, yty, SSY);
	
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
	CreateTree(branch, vars, bestmodel, model, n, m, modeldim);

	int pmodel = INTEGER(modeldim)[m];
	SEXP Rmodel_m;
	PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
	PROTECT(Rcoef_m = NEW_NUMERIC(pmodel)); 
	PROTECT(Rse_m = NEW_NUMERIC(pmodel));  

	model_m = GetModel_m(Rmodel_m, model, p);
	//evaluate logmargy and shrinkage

	R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel, p, nobs, m, &mse_m);
	gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
	
	prior_m  = compute_prior_probs(model,pmodel,p, modelprior);
	REAL(Rbestmarg)[0] = REAL(logmarg)[m];
	
	SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
	SetModel(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2, m);

	int nUnique=0, newmodel=0;
	double *real_model = vecalloc(n);
	int *modelold = ivecalloc(p);
	int old_loc = 0;
	int new_loc;
	pmodel_old = pmodel;
	nUnique=1;
	INTEGER(counts)[0] = 0;
	postold =  REAL(logmarg)[m] + log(REAL(priorprobs)[m]);
	memcpy(modelold, model, sizeof(int)*p);
	m = 0;
	int *varin= ivecalloc(p);
	int *varout= ivecalloc(p);
	double problocal = REAL(plocal)[0];
	//	while (nUnique < k && m < INTEGER(BURNIN_Iterations)[0]) {
	while (nUnique < k && m < k ) {
	        memcpy(model, modelold, sizeof(int)*p);
		pmodel =  n_sure;

		MH = GetNextModelCandidate(pmodel_old, n, n_sure, model, vars, problocal, varin, varout);
		
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
			PROTECT(Rcoef_m = NEW_NUMERIC(pmodel)); 
			PROTECT(Rse_m = NEW_NUMERIC(pmodel));   
			model_m = GetModel_m(Rmodel_m, model, p);

			R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel, p, nobs, m, &mse_m);
			gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);

			prior_m = compute_prior_probs(model,pmodel,p, modelprior);
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
				
				//record model data
				SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, nUnique);
				SetModel(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2,nUnique);
				
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
	compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);        

	//  Now sample W/O Replacement 
	//	Rprintf("NumUnique Models Accepted %d \n", nUnique);
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
	SET_STRING_ELT(ANS_names, 7, mkChar("ols"));

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
	SET_STRING_ELT(ANS_names, 8, mkChar("ols.se"));

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
	SET_STRING_ELT(ANS_names, 13, mkChar("probs.MCMC"));

	SET_VECTOR_ELT(ANS, 14, NumUnique);
	SET_STRING_ELT(ANS_names, 14, mkChar("n.Unique"));

	setAttrib(ANS, R_NamesSymbol, ANS_names);
	
	PutRNGstate();
    UNPROTECT(nProtected);
    //	Rprintf("Return\n");
	return(ANS);  
}

