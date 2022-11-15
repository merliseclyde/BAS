#include "bas.h"


/*typedef struct glmfamilystruc {
	const char *family;
	const char *link;
	void (*mu_eta)(double *eta, double *mu, int n);
	void (*linkfun)(double *mu, double *eta, int n);
	void (*variance)(double * mu, double *var, int n);
	void (*dev_resids)(double *y, double *mu, double *weights, double *resids, int n);
	void (*linkinv)(double *eta, double *mu, int n);
	void (*initialize)(double *Y, double *mu, double *weights, int n);
	double (*dispersion)(double *resid,  double *weights, int n, int rank);
} glmstptr;

added to family.h
*/
/* Version of glm.fit that can be called directly from R or C*/
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

	double  one = 1.0,  tol, devold, devnew;
	double disp;
	
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

		if (n < rank) {
			Rprintf("X has rank %ld but there are only %ld observations");
			conv = 1;
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
