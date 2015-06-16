#include "sampling.h"
#define LOG2PI  1.837877066409345
extern double hyperg1F1(double, double, double);
extern double hyperg(double, double, double);
SEXP gglm_lpy(SEXP RX, SEXP RY,SEXP Ra, SEXP Rb, SEXP Rs, SEXP Rcoef, SEXP Rmu) {
	int *xdims = INTEGER(getAttrib(RX,R_DimSymbol));
	int n=xdims[0], p = xdims[1];
	int nProtected = 0;

	SEXP ANS = PROTECT(allocVector(VECSXP, 3)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 3)); ++nProtected;
	
	//input, read only 
	double *X=REAL(RX), *Y=REAL(RY), *coef=REAL(Rcoef), *mu=REAL(Rmu);
	double a = REAL(Ra)[0], b = REAL(Rb)[0], s = REAL(Rs)[0];

	//working variables (do we really need to make them R variables?)
	SEXP RXc = PROTECT(allocVector(REALSXP,n*p)); ++nProtected;
	SEXP RIeta =  PROTECT(allocVector(REALSXP,n)); ++nProtected;  
	SEXP RXIeta=PROTECT(allocVector(REALSXP,p)); ++nProtected;
	SEXP RXcBeta =  PROTECT(allocVector(REALSXP,n)); ++nProtected;
	double *Xc=REAL(RXc), *Ieta = REAL(RIeta), *XcBeta = REAL(RXcBeta), *XIeta = REAL(RXIeta);

	//output
	SEXP RlpY=PROTECT(allocVector(REALSXP,1)); ++nProtected; 
	double lpY = NA_REAL;
	
	SEXP RQ=PROTECT(allocVector(REALSXP,1)); ++nProtected; 
	double Q = NA_REAL;

	double lC = 0.0;
	double sum_Ieta = 0.0;
	for (int i = 0; i < n; i++) {
		lC += dbinom(Y[i],1.0,mu[i],1);
		Ieta[i] = mu[i] * (1.0 - mu[i]);
		sum_Ieta += Ieta[i];
	}

	if (p == 0) { //if null model
		if ( b == 0) {
			// for Jeffreys's prior (b = 0), set lpY to NA, since then CH g-prior isn't appropriate for the null model in this case
			lpY = NA_REAL;
		} else {	
			lpY = lC + 0.5 * LOG2PI - 0.5 * log(sum_Ieta);    
		}
    } else { //not null model
		// "centering"
		for (int i = 0; i < p; i++) {
			double temp = 0.0;
			int base = i * n;
			for (int j = 0; j < n; j++) {
				temp += X[base + j] * Ieta[j];
			}
			XIeta[i] = temp / sum_Ieta;
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

		//invIbeta <- solve(t(Xc) %*% sweep(Xc,1,Ieta,FUN = "*"));
		//for (int j = 0; j < n; j++) { 
		//	Ieta[j] = sqrt(Ieta[j]); // reuse Ieta for root
		//	for (int i =0; i < p; i++) { //reuse Xc too
		//		Xc[i*n+j] *=  Ieta[j];
		//	}
		//}

		
    	lpY = lC + 0.5 * LOG2PI - 0.5 * log(sum_Ieta);
	Rprintf("log(sum_Ieta = %lf\n", log(sum_Ieta));
	lpY += lbeta((a + p) / 2.0, b / 2.0) +
	       log(hyperg1F1((a + p)/2.0, (a + b + p)/2.0, -(s + Q)/2.0)); 
    	//doesn't apply for the Jeffreys prior
		if (a > 0 && b > 0) {
		    lpY = lpY - lbeta(a / 2.0, b / 2.0) -
		          log(hyperg1F1(a/2.0, (a + b)/2.0, - s/2.0));
		}
		Rprintf("logmarg = %lf\n", lpY);
	}
	REAL(RlpY)[0] = lpY;
	REAL(RQ)[0] = Q;

	SET_VECTOR_ELT(ANS, 0, RlpY);
	SET_STRING_ELT(ANS_names, 0, mkChar("lpY"));
	SET_VECTOR_ELT(ANS, 1, RQ);
	SET_STRING_ELT(ANS_names, 1, mkChar("Q"));
	SET_VECTOR_ELT(ANS, 2, RIeta);
	SET_STRING_ELT(ANS_names, 2, mkChar("Ieta"));
	setAttrib(ANS, R_NamesSymbol, ANS_names);

	UNPROTECT(nProtected);
	return(ANS);
	//return(RlpY);
}
