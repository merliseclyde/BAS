// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include "bas.h"


/* Functions for getting MLEs and Bayes Estimates for each model */

SEXP glm_FitModel(SEXP RX, SEXP RY, SEXP Rmodel_m,  //input data
                  SEXP Roffset, SEXP Rweights, glmstptr * glmfamily, SEXP Rcontrol,
                  SEXP Rlaplace,  betapriorptr * betapriorfamily) { //parameters
  
  int nprotected = 0;
  int *model_m = INTEGER(Rmodel_m);
  int pmodel = LENGTH(Rmodel_m);
  //subset the data and call the model fitting function
  int n = INTEGER(getAttrib(RX,R_DimSymbol))[0];
  double *X = REAL(RX);
  
  
  SEXP RXmodel=PROTECT(allocMatrix(REALSXP, n , pmodel)); nprotected++;
  double *Xwork = REAL(RXmodel);
  for (int j=0; j < pmodel; j++) { //subsetting matrix
    int model_m_j = model_m[j];
    memcpy(Xwork + j * n, X + model_m_j*n, sizeof(double)*n);
  }
 
  SEXP glm_MLEs = PROTECT(glm_bas(RXmodel, RY, glmfamily, Roffset, Rweights, Rcontrol));
  nprotected++;
  
  SEXP RXmodel_noIntercept=PROTECT(allocMatrix(REALSXP, n , pmodel-1)); nprotected++;
  
  if (pmodel > 1) {
    double *Xwork_noIntercept = REAL(RXmodel_noIntercept);
    memcpy(Xwork_noIntercept, Xwork + n, sizeof(double)*n*(pmodel-1));
  }
  
  //extract mu and coef and evaluate the function
  SEXP Rmu = PROTECT(duplicate(getListElement(glm_MLEs, "mu"))); nprotected++;
  SEXP Rdeviance = PROTECT(duplicate(getListElement(glm_MLEs, "deviance"))); nprotected++;
  SEXP Rcoef = PROTECT(duplicate(getListElement(glm_MLEs, "coefficients")));nprotected++;
  
  SEXP Rlpy = PROTECT(gglm_lpy(RXmodel_noIntercept, RY, Rcoef, Rmu, Rdeviance, Rweights,
                               glmfamily, betapriorfamily,  Rlaplace));
  nprotected++;
  
  SEXP ANS = PROTECT(allocVector(VECSXP, 2)); nprotected++;
  SEXP ANS_names = PROTECT(allocVector(STRSXP, 2)); nprotected++;
  
  SET_VECTOR_ELT(ANS, 0, glm_MLEs);
  SET_VECTOR_ELT(ANS, 1, Rlpy);
  SET_STRING_ELT(ANS_names, 0, mkChar("fit"));
  SET_STRING_ELT(ANS_names, 1, mkChar("lpy"));
  
  setAttrib(ANS, R_NamesSymbol, ANS_names);
  
  UNPROTECT(nprotected);
  return(ANS);
}


SEXP gglm_lpy(SEXP RX, SEXP RY, SEXP Rcoef, SEXP Rmu, SEXP Rdeviance, SEXP Rwts, 
              glmstptr * glmfamily, betapriorptr * betapriorfamily, SEXP  Rlaplace) {
  int *xdims = INTEGER(getAttrib(RX,R_DimSymbol));
  int n=xdims[0], p = xdims[1];
  int nProtected = 0;
  
  SEXP ANS = PROTECT(allocVector(VECSXP, 5)); ++nProtected;
  SEXP ANS_names = PROTECT(allocVector(STRSXP, 5)); ++nProtected;
  
  //input, read only
  double *X=REAL(RX), *Y=REAL(RY), *coef=REAL(Rcoef), *mu=REAL(Rmu), *devb=REAL(Rdeviance), *weights=REAL(Rwts);
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
  
  double loglik_mle = 0.0, temp;
  double sum_Ieta = 0.0, logdet_Iintercept;
  int i, j, l, base;
  
  loglik_mle = glmfamily->loglik(Y, mu, weights, devb[0], n);
  glmfamily->info_matrix(Y, mu, weights, Ieta, n);
  
  for (i = 0; i < n; i++) {
    sum_Ieta += Ieta[i];
  }
  
  logdet_Iintercept = log(sum_Ieta);
  
  
  for (i = 0; i < p; i++) {
    temp = 0.0;
    base = i * n;
    for (j = 0; j < n; j++) {
      temp += X[base + j] * Ieta[j];
    }
    XIeta[i] = temp / sum_Ieta;   // Xbar in i.p. space
  }
  
  //Xc <- X - rep(1,n) %*% t((t(X) %*% Ieta)) / sum.Ieta;
  for (i =0, l =0; i < p; i++) {
    temp = XIeta[i];
    for (j = 0; j < n; j++,l++) {
      Xc[l] = X[l] - temp;
    }
  }
  
  //Q <- sum((Xc %*% beta)^2 * Ieta);
  for (j = 0; j < n; j++) { //double check if this is already zero by default
    XcBeta[j] = 0.0;
  }
  
  for (i = 0, l=0; i < p; i++) {
    double beta = coef[i+1];
    for (int j = 0; j < n; j++,l++) {
      XcBeta[j] += Xc[l] * beta;
    }
  }
  
  Q = 0.0;
  for (j = 0; j < n; j++) {
    Q += XcBeta[j] * XcBeta[j] * Ieta[j];
  }
  
  
  lpY = betapriorfamily->logmarglik_fun(betapriorfamily->hyperparams, p, Q,
                                        loglik_mle, logdet_Iintercept, laplace);
  
  shrinkage_m = betapriorfamily->shrinkage_fun(betapriorfamily->hyperparams, p, Q, laplace);
  
  intercept = coef[0];
  for (i = 0; i < p; i++) {
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


SEXP glm_bas(SEXP RX, SEXP RY, glmstptr *glmfamily, SEXP Roffset, SEXP Rweights, SEXP Rcontrol) {
	int   *xdims = INTEGER(getAttrib(RX,R_DimSymbol)), n=xdims[0], p = xdims[1];
	int inc = 1, nProtected = 0, it=0;

	SEXP ANS = PROTECT(allocVector(VECSXP, 6)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 6)); ++nProtected;
	SEXP RXwork = PROTECT(duplicate(RX)); ++nProtected;
	SEXP RYwork = PROTECT(duplicate(RY)); ++nProtected;
	SEXP RWwork = PROTECT(duplicate(RY)); ++nProtected;
	SEXP Rvariance = PROTECT(duplicate(RY)); ++nProtected;
	SEXP Rmu_eta = PROTECT(duplicate(RY)); ++nProtected;
	SEXP Reta = PROTECT(duplicate(RY)); ++nProtected;
	SEXP Rmu = PROTECT(duplicate(RY)); ++nProtected;
	SEXP Rcoef= PROTECT(allocVector(REALSXP,p)); ++nProtected;
	SEXP Rcoefwork= PROTECT(allocVector(REALSXP,p)); ++nProtected;
	SEXP Rrank=PROTECT(allocVector(INTSXP,1)); ++nProtected;
	SEXP Rcov = PROTECT(allocVector(REALSXP, p*p)); ++nProtected;
	SEXP RR = PROTECT(allocVector(REALSXP, p*p)); ++nProtected;
	SEXP Rse= PROTECT(allocVector(REALSXP, p)); ++nProtected;
	SEXP Rresiduals= PROTECT(duplicate(RY)); ++nProtected;
	SEXP Reffects= PROTECT(duplicate(RY)); ++nProtected;
	SEXP Rpivot=PROTECT(allocVector(INTSXP,p)); ++nProtected;
	SEXP Rqrauxmat=PROTECT(allocVector(REALSXP,p)); ++nProtected;
	SEXP Rworkmat=PROTECT(allocVector(REALSXP,2*p)); ++nProtected;
	SEXP Rdeviance=PROTECT(allocVector(REALSXP,1)); ++nProtected;
	SEXP RregSS=PROTECT(allocVector(REALSXP,1)); ++nProtected;

	double *X=REAL(RX), *Y=REAL(RY), *Xwork=REAL(RXwork),
		*w=REAL(RWwork),*Ywork=REAL(RYwork), *effects=REAL(Reffects),
		*coef=REAL(Rcoef),*coefwork=REAL(Rcoefwork), *se=REAL(Rse), *cov = REAL(Rcov), *R = REAL(RR),
		*work=REAL(Rworkmat), *qraux=REAL(Rqrauxmat), *weights=REAL(Rweights),
		*mu=REAL(Rmu), *offset=REAL(Roffset),*eta=REAL(Reta),  *mu_eta=REAL(Rmu_eta),
		*residuals=REAL(Rresiduals), *dev=REAL(Rdeviance), *regSS = REAL(RregSS),
		*variance=REAL(Rvariance);


	double one = 1.0,  tol, devold, devnew;
	double disp = 1.0;

	
	int   i, j, l, rank=1, *pivot=INTEGER(Rpivot), conv=0;

  // char  trans[]="N";

	//	glmstptr *glmfamily;
	//      glmfamily = make_glmfamily_structure(family);


	tol = fmin(1e-07, REAL(getListElement(Rcontrol,"epsilon"))[0]/1000);



	//fit the model
	glmfamily->initialize(Y, mu, weights, n);
	glmfamily->linkfun(mu, eta, n);
	glmfamily->linkinv(eta, mu, n);
	glmfamily->dev_resids(Y, mu, weights, residuals, n);
	devold = deviance(residuals, n);
	devnew = devold;
	conv = 0.0;
	it = 0;

	while ( conv < 1 && it < REAL(getListElement(Rcontrol, "maxit"))[0]) {
		glmfamily->mu_eta(eta, mu_eta, n);
		glmfamily->variance(mu, variance, n);

		for (i=0, l=0; i<n; i++) {
			w[i] = sqrt(weights[i]*mu_eta[i]*mu_eta[i]/variance[i]);
			Ywork[i] = w[i]*(eta[i] - offset[i] + (Y[i] - mu[i])/mu_eta[i]);
			residuals[i] = (Y[i] - mu[i])/mu_eta[i];
		}
		for (j=0, l=0; j<p; j++) {
			pivot[j] = j+1;
			for (i=0; i<n; i++, l++) {
				Xwork[l] = REAL(RX)[l]*w[i];
			}
		}

		disp = glmfamily->dispersion(residuals, weights, n, rank);
		
		
		rank = 1;
		for (j=0; j<p; j++) {
			pivot[j] = j+1;
		}

		F77_NAME(dqrls)(&Xwork[0], &n, &p, &Ywork[0], &inc, &tol,  &coefwork[0],
			&residuals[0], &effects[0], &rank, &pivot[0], &qraux[0], &work[0]);

		//    Rprintf("rank %ld \n", rank);
    //    should no get here
		if (n < rank) { // # nocov start
			error("X has rank %d but there are only %d observations", rank, n);
			conv = 1;  // # nocov end
		}

		for (j=0; j<p; j++) {
			coef[pivot[j] - 1] = coefwork[j];
		}

		F77_NAME(dcopy)(&n, &offset[0], &inc, &eta[0], &inc);
		F77_NAME(dgemv)("N", &n, &p, &one, &X[0], &n, &coef[0], &inc, &one, &eta[0],&inc 
                    FCONE);

		glmfamily->linkinv(eta, mu, n);
		glmfamily->dev_resids(Y, mu, weights, residuals, n);
		devnew = deviance(residuals, n);
		glmfamily->mu_eta(eta, mu_eta, n);
		glmfamily->variance(mu, variance, n);

		devnew = deviance(residuals, n);
		if (fabs(devnew - devold)/(0.1 + fabs(devnew)) < REAL(getListElement(Rcontrol, "epsilon"))[0]) {
			conv = 1;
		} else {
			devold = devnew;
		}
		it += 1;

		dev[0] = devnew;
		memset(se, 0.0, p*sizeof(double));
		
		if (rank == p)   {
			chol2se(&Xwork[0], &se[0], &R[0], &cov[0], p, n);
		} else {
			QR2cov(&Xwork[0], &R[0], &cov[0], rank, n);
			for (j=0; j < rank; j++) {
				se[pivot[j]-1] = cov[j*rank + j];
			}
		}
		
		for (j=0; j < p; j++) {
		  se[j] = sqrt(se[j]*disp);
		}
		

		regSS[0] = quadform(coefwork, R, rank);

		INTEGER(Rrank)[0] = rank;
	}

	SET_VECTOR_ELT(ANS, 0, Rcoef);
	SET_VECTOR_ELT(ANS, 1, Rse);
	SET_VECTOR_ELT(ANS, 2, Rmu);
	SET_VECTOR_ELT(ANS, 3, Rdeviance);
	SET_VECTOR_ELT(ANS, 4, Rrank);
	SET_VECTOR_ELT(ANS, 5, RregSS);

	SET_STRING_ELT(ANS_names, 0, mkChar("coefficients"));
	SET_STRING_ELT(ANS_names, 1, mkChar("se"));
	SET_STRING_ELT(ANS_names, 2, mkChar("mu"));
	SET_STRING_ELT(ANS_names, 3, mkChar("deviance"));
	SET_STRING_ELT(ANS_names, 4, mkChar("rank"));
	SET_STRING_ELT(ANS_names, 5, mkChar("RegSS"));

	setAttrib(ANS, R_NamesSymbol, ANS_names);

	UNPROTECT(nProtected);

	return(ANS);
}
