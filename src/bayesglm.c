#include "bas.h"

typedef struct coefpriorstruc {
  const char *family;
  const char *class;
  double *hyper;
  double (*log_marginal_likelihood)(double dev, double regSS, int n, int p, int pgamma, double g, double *hyper);
  double (*shrinkage)(double dev,  double regSS, int n, int p, int pgamma, double g, double *hyper);
  double (*g)(double dev,  double regSS, int n, int p, int pgamma, double *hyper);
} coefdistptr;




 double no_shrinkage(double dev,  double regSS, int n, int p, int pgamma, double g,  double *hyper) {
  return 1.0;
}

 double shrinkage_gprior(double dev,  double regSS, int n, int p, int pgamma,  double g, double *hyper) {
  return g/(1.0 + g) ;
}

 double g_EB_local(double dev,  double regSS, int n, int p, int pgamma, double *hyper) {
  double g = regSS/pgamma - 1.0;
  return g > 0.0 ? g : 0.0;
}

 double g_gprior(double dev,  double regSS, int n, int p, int pgamma, double *hyper) {
  return hyper[0];
}

 double no_g(double dev,  double regSS, int n, int p, int pgamma, double *hyper) {
  return 1.0;
}

double log_marginal_likelihood_IC(double dev, double regSS, int n, int p, int pgamma, double g, double *hyper) {
  return -.5*(dev +  pgamma*hyper[0]);
}
double log_marginal_likelihood_gprior(double dev, double regSS, int n, int p, int pgamma, double g, double *hyper) {
  return -.5*(dev  + pgamma*log(g + 1.0) + regSS/(g + 1.0));
}


/* Version of glm.fit that can be called directly from R or C*/

// [[register]]
SEXP glm_fit(SEXP RX, SEXP RY,SEXP family, SEXP Roffset, SEXP Rweights, SEXP Rpriorcoef, SEXP Rcontrol)
{
  int   *xdims = INTEGER(getAttrib(RX,R_DimSymbol)), n=xdims[0], p = xdims[1];
  int inc = 1,  nmodels=1,  nProtected = 0, it=0;
    // int job = 01, info = 1;


  SEXP ANS = PROTECT(allocVector(VECSXP, 9)); ++nProtected;
  SEXP ANS_names = PROTECT(allocVector(STRSXP, 9)); ++nProtected;
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
  SEXP Rlog_marg_lik=PROTECT(allocVector(REALSXP,nmodels)); ++nProtected;
  SEXP Rdeviance=PROTECT(allocVector(REALSXP,nmodels)); ++nProtected;
  SEXP RregSS=PROTECT(allocVector(REALSXP,nmodels)); ++nProtected;
  SEXP Rg = PROTECT(allocVector(REALSXP, nmodels)); ++nProtected;
  SEXP Rshrinkage = PROTECT(allocVector(REALSXP, nmodels)); ++nProtected;


  double *X=REAL(RX), *Y=REAL(RY), *Xwork=REAL(RXwork),
    *w=REAL(RWwork),*Ywork=REAL(RYwork), *effects=REAL(Reffects),
    *coef=REAL(Rcoef),*coefwork=REAL(Rcoefwork), *se=REAL(Rse), *cov = REAL(Rcov), *R = REAL(RR),
    *work=REAL(Rworkmat), *qraux=REAL(Rqrauxmat), *weights=REAL(Rweights),
    *mu=REAL(Rmu), *offset=REAL(Roffset),*eta=REAL(Reta),  *mu_eta=REAL(Rmu_eta),
    *residuals=REAL(Rresiduals), *dev=REAL(Rdeviance), *regSS = REAL(RregSS), *g = REAL(Rg),
    *variance=REAL(Rvariance), *hyper;

  double  one = 1.0,  tol, devold, devnew;

  int   i, j, l, m, rank=1, *pivot=INTEGER(Rpivot), conv=0;

  glmstptr *glmfamily;
  coefdistptr *coefprior;
  char  trans[]="N";

  tol = fmin(1e-07, REAL(getListElement(Rcontrol,"epsilon"))[0]/1000);

  coefprior = (struct coefpriorstruc *) R_alloc(1, sizeof(struct coefpriorstruc));
  coefprior->family = CHAR(STRING_ELT(getListElement(Rpriorcoef, "family"),0));
  //  Rprintf("prior family %s\n", coefprior->family);
  coefprior->class  = CHAR(STRING_ELT(getListElement(Rpriorcoef, "class"),0));
  //  Rprintf("class %s\n", coefprior->class);
  //  if (getListElement(Rpriorcoef, "hyper") != R_NilValue)
  hyper = REAL(getListElement(Rpriorcoef, "hyper"));

  if  (strcmp(coefprior->class, "gprior") == 0) {
    coefprior->shrinkage = shrinkage_gprior;
    coefprior->log_marginal_likelihood = log_marginal_likelihood_gprior;
    if  (strcmp(coefprior->family, "fixed-g-prior") == 0)
      coefprior->g = g_gprior;
    else
      coefprior->g = g_EB_local;
}
  if  (strcmp(coefprior->class, "IC") == 0) {
    coefprior->shrinkage = no_shrinkage;
    coefprior->log_marginal_likelihood = log_marginal_likelihood_IC;
    coefprior->g = no_g;
}



  glmfamily = make_glmfamily_structure(family);



  for (m=0; m< nmodels; m++){
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

    rank = 1;
    for (j=0; j<p; j++) {
      pivot[j] = j+1;
    }

    F77_NAME(dqrls)(&Xwork[0], &n, &p, &Ywork[0], &inc, &tol,  &coefwork[0],
	 &residuals[0], &effects[0], &rank, &pivot[0], &qraux[0], &work[0]);

    //    Rprintf("rank %ld \n", rank);

    if (n < rank) {
      warning("X has rank %ld but there are only %ld observations");
      conv = 1;
    }

    for (j=0; j<p; j++) {
      coef[pivot[j] - 1] = coefwork[j];
    }


    F77_NAME(dcopy)(&n, &offset[0], &inc, &eta[0], &inc);
    F77_NAME(dgemv)(trans, &n, &p, &one, &X[0], &n, &coef[0], &inc, &one, &eta[0],&inc);

    glmfamily->linkinv(eta, mu, n);
    glmfamily->dev_resids(Y, mu, weights, residuals, n);
    devnew = deviance(residuals, n);
    glmfamily->mu_eta(eta, mu_eta, n);
    glmfamily->variance(mu, variance, n);

    devnew = deviance(residuals, n);
    //    Rprintf("old %f new %f conv %d\n", devold,devnew, conv);

    if (fabs(devnew - devold)/(0.1 + fabs(devnew)) < REAL(getListElement(Rcontrol, "epsilon"))[0]) {
     conv = 1;
	  }
    else { devold=devnew;}
    it += 1;
  }

  dev[m] = devnew;


      if (rank == p)   chol2se(&Xwork[0], &se[0], &R[0], &cov[0], p, n);
      else	{
	QR2cov(&Xwork[0], &R[0], &cov[0], rank, n);
	for (j=0; j < rank; j++)  se[pivot[j]-1] = sqrt(cov[j*rank + j]);
      }

  regSS[m] = quadform(coefwork, R, rank);
  g[m] = coefprior->g(dev[m],  regSS[m],  n,  p,  rank, hyper);
  REAL(Rshrinkage)[m]  = coefprior->shrinkage(dev[m],  regSS[m],  n,  p, rank, g[m],  hyper);
  REAL(Rlog_marg_lik)[m]  = coefprior->log_marginal_likelihood(dev[m],  regSS[m],  n,  p, rank, g[m],  hyper);

  INTEGER(Rrank)[m] = rank;
  }

  SET_VECTOR_ELT(ANS, 0, Rcoef);
  SET_VECTOR_ELT(ANS, 1, Rse);
  SET_VECTOR_ELT(ANS, 2, Rmu);
  SET_VECTOR_ELT(ANS, 3, Rdeviance);
  SET_VECTOR_ELT(ANS, 4, Rrank);
  SET_VECTOR_ELT(ANS, 5, Rg);
  SET_VECTOR_ELT(ANS, 6, Rshrinkage);
  SET_VECTOR_ELT(ANS, 7, RregSS);
  SET_VECTOR_ELT(ANS, 8, Rlog_marg_lik);

  SET_STRING_ELT(ANS_names, 0, mkChar("coefficients"));
  SET_STRING_ELT(ANS_names, 1, mkChar("se"));
  SET_STRING_ELT(ANS_names, 2, mkChar("mu"));
  SET_STRING_ELT(ANS_names, 3, mkChar("deviance"));
  SET_STRING_ELT(ANS_names, 4, mkChar("rank"));
  SET_STRING_ELT(ANS_names, 5, mkChar("g"));
  SET_STRING_ELT(ANS_names, 6, mkChar("shrinkage"));
  SET_STRING_ELT(ANS_names, 7, mkChar("RegSS"));
  SET_STRING_ELT(ANS_names, 8, mkChar("logmarglik"));

  setAttrib(ANS, R_NamesSymbol, ANS_names);

  UNPROTECT(nProtected);

  return(ANS);
  }


