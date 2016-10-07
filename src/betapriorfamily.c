#include <Rinternals.h>
#include <Rconfig.h>
#include <Rmath.h>
#include <R_ext/Constants.h>
#include "sampling.h"
#include "betapriorfamily.h"
#include <float.h>

// WARNING pmodel below has intercept removed !!!  //

struct betapriorfamilystruc * make_betaprior_structure(SEXP betaprior, SEXP glmfamily) {

  betapriorptr *betapriorfamily;

  //  Rprintf("Make family\n");
  
  betapriorfamily = (struct betapriorfamilystruc *) R_alloc(1, sizeof(struct betapriorfamilystruc)) ;
  betapriorfamily->priorfamily = CHAR(STRING_ELT(getListElement(betaprior, "family"),0));
  //  Rprintf("family %s\n", betapriorfamily->priorfamily);

  betapriorfamily->priorclass = CHAR(STRING_ELT(getListElement(betaprior, "class"),0));
  //  Rprintf("family %s\n", betapriorfamily->priorclass);

  betapriorfamily->samplingmodel = CHAR(STRING_ELT(getListElement(glmfamily, "family"),0));
  // Rprintf("samplingmodel %s\n", betapriorfamily->samplingmodel);
  // SEXP hyperparameters = PROTECT(duplicate(getListElement(betaprior, "hyper.parameters")));	
  betapriorfamily->hyperparams = getListElement(betaprior, "hyper.parameters");
  
  // Rprintf("create structure \n");
  
  if (strcmp(betapriorfamily->priorfamily, "CCH") == 0) {
    betapriorfamily->logmarglik_fun = CCH_glm_logmarg;
    betapriorfamily->shrinkage_fun = CCH_glm_shrinkage;
  }
  else if (strcmp(betapriorfamily->priorfamily, "tCCH") == 0) {
    betapriorfamily->logmarglik_fun = tCCH_glm_logmarg;
    betapriorfamily->shrinkage_fun = tCCH_glm_shrinkage;  
  }
    else if (strcmp(betapriorfamily->priorfamily, "intrinsic") == 0) {
    betapriorfamily->logmarglik_fun = intrinsic_glm_logmarg;
    betapriorfamily->shrinkage_fun = intrinsic_glm_shrinkage;  
  }

  else if (strcmp(betapriorfamily->priorfamily, "hyper-g/n") == 0) {
    betapriorfamily->logmarglik_fun = tCCH_glm_logmarg;
    betapriorfamily->shrinkage_fun = tCCH_glm_shrinkage;  
  }
  else if (strcmp(betapriorfamily->priorfamily, "Jeffreys") == 0) {
    betapriorfamily->logmarglik_fun = Jeffreys_glm_logmarg;
    betapriorfamily->shrinkage_fun = CCH_glm_shrinkage;  // OK
  }

  else  if (strcmp(betapriorfamily->priorclass, "IC") == 0) {
    betapriorfamily->logmarglik_fun = IC_glm_logmarg;
    betapriorfamily->shrinkage_fun = IC_shrinkage;
  }
  else if (strcmp(betapriorfamily->priorfamily, "robust") == 0) {
    betapriorfamily->logmarglik_fun = robust_glm_logmarg;
    betapriorfamily->shrinkage_fun = robust_glm_shrinkage;
    //    Rprintf("n = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "n"))[0]);
  }

  else if (strcmp(betapriorfamily->priorfamily, "betaprime") == 0) {
     betapriorfamily->logmarglik_fun = betaprime_glm_logmarg;
     betapriorfamily->shrinkage_fun = betaprime_glm_shrinkage;
     //     Rprintf("n = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "n"))[0]);
     //    Rprintf("a = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "alpha"))[0]);
  }

  else if (strcmp(betapriorfamily->priorfamily, "TG") == 0) {
     betapriorfamily->logmarglik_fun = TG_glm_logmarg;
     betapriorfamily->shrinkage_fun = TG_glm_shrinkage;
     //     Rprintf("a = %lf\n", REAL(getListElement(betapriorfamily->hyperparams, "alpha"))[0]);
  }

  else if (strcmp(betapriorfamily->priorfamily, "EB-local") == 0) {
    betapriorfamily->logmarglik_fun = EB_local_glm_logmarg;
    betapriorfamily->shrinkage_fun = EB_local_glm_shrinkage;
  }

  else if (strcmp(betapriorfamily->priorfamily, "g.prior") == 0) {
    betapriorfamily->logmarglik_fun = g_prior_glm_logmarg;
    betapriorfamily->shrinkage_fun = g_prior_shrinkage;
  }
  else if (strcmp(betapriorfamily->priorfamily, "testBF.prior") == 0) {
    betapriorfamily->logmarglik_fun = testBF_prior_glm_logmarg;
    betapriorfamily->shrinkage_fun = g_prior_shrinkage;
  }
  else error("Prior %s has not been implemented or is misspelled\n", betapriorfamily->priorfamily);
  return(betapriorfamily);
}



double CCH_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double a, b, s, logmarglik, p;
   
  a = REAL(getListElement(hyperparams, "alpha"))[0];
  b = REAL(getListElement(hyperparams, "beta"))[0];
  s = REAL(getListElement(hyperparams, "s"))[0];
  //  n = INTEGER(getListElement(hyperparams, "n"))[0];
  //  p = INTEGER(getListElement(hyperparams, "p"))[0];
  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  p = (double) pmodel;

  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (p >= 1.0) {
    logmarglik +=   lbeta((a + p) / 2.0, b / 2.0) 
                  + loghyperg1F1((a + p)/2.0, (a + b + p)/2.0, -(s+W)/2.0, Laplace)
                  - lbeta(a / 2.0, b / 2.0)
                  - loghyperg1F1(a/2.0, (a + b)/2.0, - s/2.0, Laplace);
  }

  return(logmarglik);
}
double Jeffreys_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double a, b, s, logmarglik, p;
   
  a = REAL(getListElement(hyperparams, "alpha"))[0];
  b = REAL(getListElement(hyperparams, "beta"))[0];
  s = REAL(getListElement(hyperparams, "s"))[0];
  //  n = INTEGER(getListElement(hyperparams, "n"))[0];
  //  p = INTEGER(getListElement(hyperparams, "p"))[0];
  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  p = (double) pmodel;

  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (p >= 1.0) {
    logmarglik +=   lbeta((a + p) / 2.0, b / 2.0) 
      + loghyperg1F1((a + p)/2.0, (a + b + p)/2.0, -(s+W)/2.0, Laplace);
  }

  return(logmarglik);
}


double tCCH_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double a, b, s, r, v, theta, logmarglik, p;
   
  a = REAL(getListElement(hyperparams, "alpha"))[0];
  b = REAL(getListElement(hyperparams, "beta"))[0];
  s = REAL(getListElement(hyperparams, "s"))[0];
  r = REAL(getListElement(hyperparams, "r"))[0];
  v = REAL(getListElement(hyperparams, "v"))[0];
  theta =  REAL(getListElement(hyperparams, "theta"))[0];

  p = (double) pmodel;

  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (p >= 1.0) {
    logmarglik +=   lbeta((a + p) / 2.0, b / 2.0) 
      + log(HyperTwo(b/2.0, r, (a + b + p)/2.0, (s+W)/(2.0*v), 1.0 - theta))
      -.5*p*log(v) -.5*W/v
      - lbeta(a / 2.0, b / 2.0)
      - log(HyperTwo(b/2.0, r, (a + b)/2.0, s/(2.0*v), 1.0 - theta));
  }	

  return(logmarglik);
}

double tCCH_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double a, b, s, r, v, theta, p, shrinkage;
  
  a = REAL(getListElement(hyperparams, "alpha"))[0];
  b = REAL(getListElement(hyperparams, "beta"))[0];
  s = REAL(getListElement(hyperparams, "s"))[0];
  r = REAL(getListElement(hyperparams, "r"))[0];
  v = REAL(getListElement(hyperparams, "v"))[0];
  theta =  REAL(getListElement(hyperparams, "theta"))[0];
  
  p = (double) pmodel;
  
  shrinkage = 1.0;
  if (p >= 1.0) {
   shrinkage -=  exp( -log(v)
    + lbeta((a + p) / 2.0 + 1.0, b / 2.0) 
    + log(HyperTwo(b/2.0, r, (a +b+p)/2.0 + 1.0, (s+W)/(2.0*v), 1.0-theta)) 
    - lbeta((a+p) / 2.0, b/2.0)
    - log(HyperTwo(b/2.0, r, (a + p+ b)/2.0, (s+W)/(2.0*v), 1.0 - theta)));
  }	
  
  return(shrinkage);
}

double intrinsic_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double a, b, s, r, v, theta,n, logmarglik, p;
   
  a = REAL(getListElement(hyperparams, "alpha"))[0];
  b = REAL(getListElement(hyperparams, "beta"))[0];
  s = REAL(getListElement(hyperparams, "s"))[0];
  r = REAL(getListElement(hyperparams, "r"))[0];
  n = REAL(getListElement(hyperparams, "n"))[0];

  p = (double) pmodel;

  v = (n + p + 1.0)/(p + 1);
  theta = (n + p + 1.0)/n;
  
  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (p >= 1.0) {
    logmarglik +=   lbeta((a + p) / 2.0, b / 2.0) 
      + log(HyperTwo(b/2.0, r, (a + b + p)/2.0, (s+W)/(2.0*v), 1.0 - theta))
      -.5*p*log(v) -.5*W/v
      - lbeta(a / 2.0, b / 2.0)
      - log(HyperTwo(b/2.0, r, (a + b)/2.0, s/(2.0*v), 1.0 - theta));
  }	

  return(logmarglik);
}

double intrinsic_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double a, b, s, r, v, theta, n, p, u, shrinkage;
  
  a = REAL(getListElement(hyperparams, "alpha"))[0];
  b = REAL(getListElement(hyperparams, "beta"))[0];
  s = REAL(getListElement(hyperparams, "s"))[0];
  r = REAL(getListElement(hyperparams, "r"))[0];
  n = REAL(getListElement(hyperparams, "n"))[0];
    
  p = (double) pmodel;
  v = (n + p + 1.0)/(p + 1);
  theta = (n + p + 1.0)/n;
  
  shrinkage = 1.0;
  if (p >= 1.0) {
     u = exp(-log(v)
             + lbeta((a + p) / 2.0 + 1.0, b / 2.0) 
             + log(HyperTwo(b/2.0, r, (a +b+p)/2.0 + 1.0, (s+W)/(2.0*v), 1.0-theta)) 
             - lbeta((a+p) / 2.0, b/2.0)
             - log(HyperTwo(b/2.0, r, (a + p+ b)/2.0, (s+W)/(2.0*v), 1.0 - theta)));
    shrinkage = 1.0 - u;
  }	
  
  return(shrinkage);
}


double CCH_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double a, b, s, p, shrinkage = 1.0;
   
  a = REAL(getListElement(hyperparams, "alpha"))[0];
  b = REAL(getListElement(hyperparams, "beta"))[0];
  s = REAL(getListElement(hyperparams, "s"))[0];

  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  // Rprintf("s = %lf\n", s);

  p = (double) pmodel;
  shrinkage = 1.0;
  if (p >= 1.0) 
    // shrinkage = shrinkage_chg(a + p, a + b + p, -(s+W), Laplace);
    shrinkage = 1.0 - exp(log(a + p) -log(a + b + p) 
                    +loghyperg1F1((a+p+2.0)/2.0, (a+p+b +2.0)/2.0, -(s+W)/2, Laplace)
                    -loghyperg1F1((a+p)/2,(a+p+b)/2.0, -(s+W)/2, Laplace)
    );
      
  return(shrinkage);
}

double betaprime_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double a, n, p, logmarglik;

  a = REAL(getListElement(hyperparams, "alpha"))[0];
  n = REAL(getListElement(hyperparams, "n"))[0];
  p = (double) pmodel;
  
  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (p >= 1.0) {
    logmarglik +=   lbeta((a + p) / 2.0, (n - p - 1.5) / 2.0) 
      + loghyperg1F1((a + p)/2.0, (a + n - 1.5)/2.0, -W/2.0, Laplace)
      - lbeta(a / 2.0, (n - p - 1.5)/ 2.0)
      - loghyperg1F1(a/2.0, (a + n - p - 1.5)/2.0, 0.0, Laplace);
  }
  
  return(logmarglik);
}

double betaprime_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double a, n,p, b, shrinkage = 1.0;
   

  a = REAL(getListElement(hyperparams, "alpha"))[0];
  n = REAL(getListElement(hyperparams, "n"))[0];
  p = (double) pmodel;
  b = n - p - 1.5;
  
  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  // Rprintf("s = %lf\n", s);
  shrinkage = 1.0;
  if (p >= 1.0) 
    if (p >= 1.0) 
      // shrinkage = shrinkage_chg(a + p, a + b + p, -(s+W), Laplace);
      shrinkage = 1.0 - exp(log(a + p) -log(a + b + p) 
                            +loghyperg1F1((a+p+2.0)/2.0, (a+p+b +2.0)/2.0, -W/2.0, Laplace)
                            -loghyperg1F1((a+p)/2,(a+p+b)/2.0, -W/2.0, Laplace)
      );
  return(shrinkage);
}


double robust_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double n, p, logmarglik;
   
  n = REAL(getListElement(hyperparams, "n"))[0];
  p = (double) pmodel;
  
  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (p >= 1.0) {
    logmarglik += -log(2.0) + 0.5 *(log(n + 1.0) - log(p + 1.0)) 
                  +  lgammafn((p+1.0)/2.0)
                  - .5*(p + 1.0)*log(W/2.0) +
                  pgamma((p + 1.0)/(n + 1.0), 0.5*(p+1.0), 2.0/W, 1, 1);
  }
  return(logmarglik);
}


double robust_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double n, p, shrinkage = 1.0;
   
  n = REAL(getListElement(hyperparams, "n"))[0];

  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  // Rprintf("s = %lf\n", s);
  p = (double) pmodel;
  
  if (pmodel >= 1.0) {
    shrinkage = 1.0 - exp(log(0.5*(p+1.0)) - log(W/2.0) +
			  pgamma((p+1.0)/(n + 1.0), 0.5*(p+1.0) + 1.0, 2.0/W, 1, 1) -
			  pgamma((p+1.0)/(n + 1.0), 0.5*(p+1.0), 2.0/W, 1, 1));
  }
  return(shrinkage);
}

double TG_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double a,p, logmarglik;

  a = REAL(getListElement(hyperparams, "alpha"))[0];
  p = (double) pmodel;
  
  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (pmodel >= 1.0) {
    logmarglik +=   -log(2.0) + log(a)
      + 	lgammafn((a + p)/2.0)
      - 	.5*(a + p)*log(W/2.0)
      + 	pgamma(1.0, .5*(a + p), 2.0/W, 1, 1);
  }
  return(logmarglik);
}


double TG_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double a, p, shrinkage = 1.0;


  a = REAL(getListElement(hyperparams, "alpha"))[0];
  p = (double) pmodel;
  
  shrinkage = 1.0;
  if (p >= 1.0) {
    shrinkage = 1.0 - exp(log(a + p) - log(.5*W) +
			  pgamma(1.0, .5*(a + p) + 1.0, 2.0/W, 1, 1) - 
			  pgamma(1.0, .5*(a + p), 2.0/W, 1, 1));
  }
  return(shrinkage);
}

double EB_local_glm_logmarg(SEXP hyperparams, int pmodel, double W,
			    double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double ghat, p, logmarglik;
  p = (double) pmodel;

  ghat = fmax(0.0, W/p - 1);
  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (p >= 1.0 & ghat > 0) {
    logmarglik +=   -.5*p*log(1.0 + ghat) -.5*W/(1.0 + ghat);
  }
  return(logmarglik);
}


double EB_local_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  double ghat,p, shrinkage = 1.0;

  p = (double) pmodel;
  //  use_deviance = REAL(getListElement(hyperparams, "use_deviance"))[0];

  // Rprintf("a = %lf\n", a);
  // Rprintf("b = %lf\n", b);
  // Rprintf("s = %lf\n", s);

  if (p >= 1.0) {
    ghat = fmax(0.0, W/p - 1);
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

double testBF_prior_glm_logmarg(SEXP hyperparams, int pmodel, double W,
                             double loglik_mle, double logdet_Iintercept,
                             int Laplace ) {
    double g, logmarglik, loglik_null, z;
    
    g = REAL(getListElement(hyperparams, "g"))[0];
        
    loglik_null = REAL(getListElement(hyperparams, "loglik_null"))[0];

    // pmodel is 0 for null model      
    z = -2.0*(loglik_mle - loglik_null);
    logmarglik = - 0.5* (((double) pmodel)*(log(1.0 + g))  + z*g/(1.0 + g));
//    Rprintf("z = %lf  p = %d \n", z, pmodel );
    return(logmarglik);
  } 
  
double g_prior_glm_logmarg(SEXP hyperparams, int pmodel, double W,
		       double loglik_mle, double logdet_Iintercept, int Laplace ) {
  double g, logmarglik;
   
   g = REAL(getListElement(hyperparams, "g"))[0];

  logmarglik =   loglik_mle + M_LN_SQRT_2PI - 0.5* logdet_Iintercept;
  if (pmodel >= 1.0) {
    logmarglik +=   -.5*((double) pmodel)*log(1.0 + g) -.5*W/(1.0 + g);
  }
  return(logmarglik);
}


double g_prior_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace ) {
  
  double g, shrinkage;
  g = REAL(getListElement(hyperparams, "g"))[0];
  if (pmodel >= 1) {
    shrinkage = g/(1.0 + g);}
  else {
    shrinkage = 0.0;}
  
  return(shrinkage);
}
