#include "sampling.h"
#include "family.h"
#include "bas-glm.h"


	
SEXP glm_FitModel(SEXP RX, SEXP RY, SEXP Rmodel_m,  //input data
			  SEXP Roffset, SEXP Rweights, glmstptr * glmfamily, SEXP Rcontrol,
			  SEXP Ra, SEXP Rb, SEXP Rs) { //parameters
  int nprotected = 0;
  int *model_m = INTEGER(Rmodel_m);
  int pmodel = LENGTH(Rmodel_m);
	//subset the data and call the model fitting function
  int n = INTEGER(getAttrib(RX,R_DimSymbol))[0];
  double *X = REAL(RX);

  
  SEXP RXnow=PROTECT(allocMatrix(REALSXP, n , pmodel)); nprotected++;
  double *Xwork = REAL(RXnow);
  for (int j=0; j < pmodel; j++) { //subsetting matrix
      int model_m_j = model_m[j];
      memcpy(Xwork + j * n, X + model_m_j*n, sizeof(double)*n);
    }
  SEXP glm_fit = PROTECT(glm_bas(RXnow, RY, glmfamily, Roffset, Rweights, Rcontrol));
  nprotected++;
	
    //extract mu and coef and evaluate the function
  SEXP Rmu = PROTECT(duplicate(getListElement(glm_fit, "mu"))); nprotected++;
  SEXP Rcoef = PROTECT(duplicate(getListElement(glm_fit, "coefficients")));nprotected++;
  SEXP RXnow_noIntercept=PROTECT(allocMatrix(REALSXP, n , pmodel-1)); nprotected++;
  if (pmodel > 1) {
    double *Xwork_noIntercept = REAL(RXnow_noIntercept);
    memcpy(Xwork_noIntercept, Xwork + n, sizeof(double)*n*(pmodel-1));
  }

  
  SEXP Rlpy = PROTECT(gglm_lpy(RXnow_noIntercept, RY, Ra, Rb, Rs, Rcoef, Rmu, glmfamily));
  nprotected++;
	
  SEXP ANS = PROTECT(allocVector(VECSXP, 2)); nprotected++;
  SEXP ANS_names = PROTECT(allocVector(STRSXP, 2)); nprotected++;
	
  SET_VECTOR_ELT(ANS, 0, glm_fit);
  SET_VECTOR_ELT(ANS, 1, Rlpy);
  SET_STRING_ELT(ANS_names, 0, mkChar("fit"));
  SET_STRING_ELT(ANS_names, 1, mkChar("lpy"));

  setAttrib(ANS, R_NamesSymbol, ANS_names);

  UNPROTECT(nprotected);
  return(ANS);
}

SEXP glm_mcmc(SEXP Y, SEXP X, SEXP Roffset, SEXP Rweights, 
			  SEXP Rprobinit, SEXP Rmodeldim, 
			  SEXP modelprior, SEXP Rbestmodel,  SEXP plocal, 
			  SEXP BURNIN_Iterations,
			  SEXP Ra, SEXP Rb, SEXP Rs,
			  SEXP family, SEXP Rcontrol
			  )
{
	int nProtected = 0;
	int nModels=LENGTH(Rmodeldim);
	SEXP ANS = PROTECT(allocVector(VECSXP, 16)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 16)); ++nProtected;
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
	
	SEXP NumUnique = PROTECT(allocVector(INTSXP, 1)); ++nProtected;
	
	double *probs, MH=0.0, prior_m=1.0, shrinkage_m, logmargy, postold, postnew;
	int i, m, n, pmodel_old, *bestmodel;
	int mcurrent, n_sure;
	glmstptr *glmfamily;

	glmfamily = make_glmfamily_structure(family);

	//get dimsensions of all variables 
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int k = LENGTH(modelprobs);
	
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
	SEXP Rmodel_m =	PROTECT(allocVector(INTSXP,pmodel));
	GetModel_m(Rmodel_m, model, p);
	//evaluate logmargy and shrinkage
	SEXP glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights,
					    glmfamily, Rcontrol, Ra, Rb, Rs));	
	prior_m  = compute_prior_probs(model,pmodel,p, modelprior);
	
	logmargy = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
	shrinkage_m = REAL(getListElement(getListElement(glm_fit, "lpy"),
					"shrinkage"))[0];
	SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
	SetModel1(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q, m);
	UNPROTECT(2);

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
	while (nUnique < k && m < INTEGER(BURNIN_Iterations)[0]) {
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
		  GetModel_m(Rmodel_m, model, p);
		  
		  glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights,
						 glmfamily, Rcontrol, Ra, Rb, Rs));	
		  prior_m = compute_prior_probs(model,pmodel,p, modelprior);
			
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
				new_loc = nUnique;
				insert_model_tree(tree, vars, n, model, nUnique);
				INTEGER(modeldim)[nUnique] = pmodel;
				//Rprintf("model %d: %d variables\n", m, pmodel);
				SetModel2(logmargy, NA_REAL, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, nUnique);
				SetModel1(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q,nUnique);
				UNPROTECT(2);	
				++nUnique; 
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
		SEXP devianceP = PROTECT(allocVector(REALSXP, nUnique));
		for (i =0; i < nUnique; i++) { 
			REAL(devianceP)[i] = REAL(deviance)[i];
		}
		SET_VECTOR_ELT(ANS, 6, devianceP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 6, deviance);
	}
	SET_STRING_ELT(ANS_names, 6, mkChar("deviance"));

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
	SET_STRING_ELT(ANS_names, 7, mkChar("coefficients"));

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
	SET_STRING_ELT(ANS_names, 8, mkChar("se"));

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

	if (nUnique < nModels) {
		SEXP QP = PROTECT(allocVector(REALSXP, nUnique));
		for (i =0; i < nUnique; i++) { 
			REAL(QP)[i] = REAL(Q)[i];
		}
		SET_VECTOR_ELT(ANS, 15, QP);
		UNPROTECT(1);
	} else {
		SET_VECTOR_ELT(ANS, 15, Q);
	}
	SET_STRING_ELT(ANS_names, 15, mkChar("Q"));

	setAttrib(ANS, R_NamesSymbol, ANS_names);
	
	PutRNGstate();
    UNPROTECT(nProtected);
	//Rprintf("Return\n");
	return(ANS);  
}

