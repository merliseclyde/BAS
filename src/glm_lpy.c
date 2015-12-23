#include "sampling.h"
#include "family.h"
#include "betapriorfamily.h"
#include "bas-glm.h"


extern double loghyperg1F1(double, double, double, int);
extern double hyperg(double, double, double);
extern double shrinkage_chg(double a, double b, double Q, int laplace);

SEXP glm_FitModel(SEXP RX, SEXP RY, SEXP Rmodel_m,  //input data
		  SEXP Roffset, SEXP Rweights, glmstptr * glmfamily, SEXP Rcontrol,
		  SEXP Rlaplace,  betapriorptr * betapriorfamily) { //parameters
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

  
  SEXP Rlpy = PROTECT(gglm_lpy(RXnow_noIntercept, RY, Rcoef, Rmu,
			       glmfamily, betapriorfamily,  Rlaplace));
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


SEXP gglm_lpy(SEXP RX, SEXP RY, SEXP Rcoef, SEXP Rmu, glmstptr * glmfamily, betapriorptr * betapriorfamily, SEXP  Rlaplace) {
	int *xdims = INTEGER(getAttrib(RX,R_DimSymbol));
	int n=xdims[0], p = xdims[1];
	int nProtected = 0;  

	SEXP ANS = PROTECT(allocVector(VECSXP, 5)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 5)); ++nProtected;
	
	//input, read only 
	double *X=REAL(RX), *Y=REAL(RY), *coef=REAL(Rcoef), *mu=REAL(Rmu);
	int laplace = INTEGER(Rlaplace)[0];
	
	//working variables (do we really need to make them R variables?)
	SEXP RXc = PROTECT(allocVector(REALSXP,n*p)); ++nProtected;
	SEXP RIeta =  PROTECT(allocVector(REALSXP,n)); ++nProtected;  
	SEXP RXIeta=PROTECT(allocVector(REALSXP,p)); ++nProtected;
	SEXP RXcBeta =  PROTECT(allocVector(REALSXP,n)); ++nProtected;
	double *Xc=REAL(RXc), *Ieta = REAL(RIeta), *XcBeta = REAL(RXcBeta), *XIeta = REAL(RXIeta);

	//output
	SEXP Rintercept=PROTECT(allocVector(REALSXP,1)); ++nProtected;
	double intercept=NA_REAL;
	
	SEXP RlpY=PROTECT(allocVector(REALSXP,1)); ++nProtected; 
	double lpY = NA_REAL;
	
	SEXP RQ=PROTECT(allocVector(REALSXP,1)); ++nProtected; 
	double Q = NA_REAL;

	SEXP Rshrinkage=PROTECT(allocVector(REALSXP,1)); ++nProtected; 
	double shrinkage_m = 1.0;

	double loglik_mle = 0.0, temp = 0.0;
	double sum_Ieta = 0.0, logdet_Iintercept;
	int i, j,l, base;
	

	loglik_mle = glmfamily->loglik(Y, mu, n);
	glmfamily->info_matrix(Y, mu, Ieta, n);

	for (int i = 0; i < n; i++) {
	        sum_Ieta += Ieta[i];
	}

	logdet_Iintercept = log(sum_Ieta);
	

	for (int i = 0; i < p; i++) {
	  double temp = 0.0;
	  int base = i * n;
	 for (int j = 0; j < n; j++) {
	   temp += X[base + j] * Ieta[j];
	 }
	 XIeta[i] = temp / sum_Ieta;   // Xbar in i.p. space
	}
       
       //Xc <- X - rep(1,n) %*% t((t(X) %*% Ieta)) / sum.Ieta;
	for (int i =0, l =0; i < p; i++) {
	   double temp = XIeta[i];
	 for (int j = 0; j < n; j++,l++) {
	   Xc[l] = X[l] - temp;
	 }
       }

       //Q <- sum((Xc %*% beta)^2 * Ieta);
       for (int j = 0; j < n; j++) { //double check if this is already zero by default
	 XcBeta[j] = 0.0;
       }
       
       for (int i = 0,l=0; i < p; i++) {
	 double beta = coef[i+1];
	 for (int j = 0; j < n; j++,l++) {
	   XcBeta[j] += Xc[l] * beta;
	 }
       }

       Q = 0.0;
       for (int j = 0; j < n; j++) { 
	 Q += XcBeta[j] * XcBeta[j] * Ieta[j];
       }

       
       lpY = betapriorfamily->logmarglik_fun(betapriorfamily->hyperparams, p, Q,
					     loglik_mle, logdet_Iintercept, laplace);

     shrinkage_m = betapriorfamily->shrinkage_fun(betapriorfamily->hyperparams, p, Q, laplace);

     intercept = coef[0];
     for ( int i = 0; i < p; i++) {
       intercept += XIeta[i]*coef[i+1]*(1.0 - shrinkage_m);
     }
     REAL(Rintercept)[0] = intercept;
     REAL(RlpY)[0] = lpY;
     REAL(RQ)[0] = Q;
     REAL(Rshrinkage)[0] = shrinkage_m;
    
 
     SET_VECTOR_ELT(ANS, 0, RlpY);
     SET_STRING_ELT(ANS_names, 0, mkChar("lpY"));
     SET_VECTOR_ELT(ANS, 1, RQ);
     SET_STRING_ELT(ANS_names, 1, mkChar("Q"));
     SET_VECTOR_ELT(ANS, 2, RIeta);
     SET_STRING_ELT(ANS_names, 2, mkChar("Ieta"));
     SET_VECTOR_ELT(ANS, 3, Rshrinkage);
     SET_STRING_ELT(ANS_names, 3, mkChar("shrinkage"));
  
     SET_VECTOR_ELT(ANS, 4, Rintercept);
     SET_STRING_ELT(ANS_names, 4, mkChar("intercept"));
     
     setAttrib(ANS, R_NamesSymbol, ANS_names);
     
     UNPROTECT(nProtected);
     return(ANS);
	//return(RlpY);
}
