#include <Rinternals.h>
#include <Rconfig.h>
#include <Rmath.h>
#include <R_ext/Constants.h>
#include "sampling.h"
#include "betapriorfamily.h"
#include <float.h>


struct betapriorfamilystruc * make_betaprior_structure(SEXP betaprior, SEXP glmfamily) {

  betapriorptr *betapriorfamily;

  //  Rprintf("Make family\n");
  
  betapriorfamily = (struct betapriorfamilystruc *) R_alloc(1, sizeof(struct betapriorfamilystruc)) ;
  betapriorfamily->priorfamily = CHAR(STRING_ELT(getListElement(betaprior, "family"),0));
  Rprintf("family %s\n", betapriorfamily->priorfamily);

  betapriorfamily->priorclass = CHAR(STRING_ELT(getListElement(betaprior, "class"),0));
  Rprintf("family %s\n", betapriorfamily->priorclass);

  betapriorfamily->samplingmodel = CHAR(STRING_ELT(getListElement(glmfamily, "family"),0));
  Rprintf("samplingmodel %s\n", betapriorfamily->samplingmodel);
  // SEXP hyperparameters = PROTECT(duplicate(getListElement(betaprior, "hyper.parameters")));	
  betapriorfamily->hyperparams = getListElement(betaprior, "hyper.parameters");
  
  Rprintf("create structure \n");
  
  if (strcmp(betapriorfamily->priorfamily, "CCH") == 0) {
    betapriorfamily->logmarglik_fun = CCH_glm_logmarg;
    betapriorfamily->shrinkage_fun = CCH_glm_shrinkage;
    Rprintf("alpha = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "alpha"))[0]);
    Rprintf("beta = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "beta"))[0]);
    Rprintf("s = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "s"))[0]);

  }

  if (strcmp(betapriorfamily->priorclass, "IC") == 0) {
    betapriorfamily->logmarglik_fun = IC_glm_logmarg;
    betapriorfamily->shrinkage_fun = IC_shrinkage;
    Rprintf("penalty = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "penalty"))[0]);
  }


   if (strcmp(betapriorfamily->priorfamily, "robust") == 0) {
    betapriorfamily->logmarglik_fun = robust_glm_logmarg;
    betapriorfamily->shrinkage_fun = robust_glm_shrinkage;
    Rprintf("n = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "n"))[0]);
  }

   if (strcmp(betapriorfamily->priorfamily, "betaprime") == 0) {
     betapriorfamily->logmarglik_fun = betaprime_glm_logmarg;
     betapriorfamily->shrinkage_fun = betaprime_glm_shrinkage;
     Rprintf("n = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "n"))[0]);
     Rprintf("a = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "alpha"))[0]);
  }

   if (strcmp(betapriorfamily->priorfamily, "EB-local") == 0) {
    betapriorfamily->logmarglik_fun = EB_local_glm_logmarg;
    betapriorfamily->shrinkage_fun = EB_local_glm_shrinkage;
  }

  return(betapriorfamily);
}


double CCH_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double a, b, s, logmarglik;
   
  a = REAL(getListElement(hyperparams, "alpha"))[0];
  b = REAL(getListElement(hyperparams, "beta"))[0];
  s = REAL(getListElement(hyperparams, "s"))[0];
  //  n = INTEGER(getListElement(hyperparams, "n"))[0];
  //  p = INTEGER(getListElement(hyperparams, "p"))[0];
  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  // Rprintf("s = %lf\n", s);
  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (pmodel > 1.0) {
    logmarglik +=   lbeta((a + pmodel -1.0) / 2.0, b / 2.0) 
                  + loghyperg1F1((a + pmodel -1.0)/2.0, (a + b + pmodel - 1.0)/2.0, -(s+W)/2.0, Laplace)
                  - lbeta(a / 2.0, b / 2.0)
                  - loghyperg1F1(a/2.0, (a + b)/2.0, - s/2.0, Laplace);
  }
  return(logmarglik);
}


double CCH_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double a, b, s, shrinkage = 1.0;
   
  a = REAL(getListElement(hyperparams, "alpha"))[0];
  b = REAL(getListElement(hyperparams, "beta"))[0];
  s = REAL(getListElement(hyperparams, "s"))[0];

  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  // Rprintf("s = %lf\n", s);

  if (pmodel > 1.0) shrinkage = shrinkage_chg(a + pmodel -1.0, a + b + pmodel -1.0, -(s+W), Laplace);
  return(shrinkage);
}

double betaprime_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double a, n, logmarglik;

  a = REAL(getListElement(hyperparams, "alpha"))[0];
  n = REAL(getListElement(hyperparams, "n"))[0];

  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (pmodel > 1.0) {
    logmarglik +=   lbeta((a + pmodel - 1.0) / 2.0, (n - (pmodel- 1.0) - 1.5) / 2.0) 
                  + loghyperg1F1((a + pmodel - 1.0 )/2.0, (a + n - 1.5)/2.0, -W/2.0, Laplace)
                  - lbeta(a / 2.0, (n - pmodel -1.0 - 1.5)/ 2.0)
                  - loghyperg1F1(a/2.0, (a + n - pmodel -1.0 - 1.5)/2.0, 0.0, Laplace);
  }
  return(logmarglik);
}

double betaprime_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double a, n, shrinkage = 1.0;
   

  a = REAL(getListElement(hyperparams, "alpha"))[0];
  n = REAL(getListElement(hyperparams, "n"))[0];

  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  // Rprintf("s = %lf\n", s);

  if (pmodel > 1.0) shrinkage = shrinkage_chg(a + pmodel -1.0, a + n - 1.5 , -W, Laplace);
  return(shrinkage);
}


double robust_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double n, logmarglik;
   
  n = REAL(getListElement(hyperparams, "n"))[0];

  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (pmodel > 1.0) {
    logmarglik +=   -log(2.0) +.5 *log(n + 1.0) - .5+log(pmodel)
      - .5*pmodel*log(W/2.0) + pgamma(pmodel/(n + 1.0), .5*(pmodel), 2.0*(n+1.0)/(W*pmodel), 1, 1);
  }
  return(logmarglik);
}


double robust_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double n, shrinkage = 1.0;
   
  n = REAL(getListElement(hyperparams, "n"))[0];

  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  // Rprintf("s = %lf\n", s);

  if (pmodel > 1.0) {
    shrinkage = 1.0 - exp(pgamma((pmodel -1.0)/(n + 1.0), .5*(pmodel -1.0) + 1.0, 2.0/W, 1, 1) -
			  pgamma((pmodel -1.0)/(n + 1.0), .5*(pmodel -1.0), 2.0/W, 1, 1));
  }
  return(shrinkage);
}

double EB_local_glm_logmarg(SEXP hyperparams, int pmodel, double W,
			    double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double ghat, logmarglik;


  ghat = fmax(0.0, W/pmodel - 1);
  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (pmodel > 1.0 & ghat > 0) {
    logmarglik +=   -.5*(pmodel -1.0)*log(1.0 + ghat) -.5*W/(1.0 + ghat);
  }
  return(logmarglik);
}


double EB_local_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double use_deviance, ghat, shrinkage = 1.0;
   
  //  use_deviance = REAL(getListElement(hyperparams, "use_deviance"))[0];

  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  // Rprintf("s = %lf\n", s);

  if (pmodel > 1.0) {
    ghat = fmax(0.0, W/pmodel - 1);
    if (ghat > 0) {
      shrinkage = ghat/(1.0 + ghat);
    }
    else {
      shrinkage = 0.0;
    }
  }
  return(shrinkage);
}


double IC_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double penalty, logmarglik;
   
  penalty = REAL(getListElement(hyperparams, "penalty"))[0];
  logmarglik =   loglik_mle - .5*penalty*pmodel;
  return(logmarglik);
}


double IC_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double shrinkage = 1.0;
  return(shrinkage);
}
