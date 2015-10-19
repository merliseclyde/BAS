/* version  5/20/2005 */
/* Rsample.c progrma for sampling without replacement in R  MC 11/2002 */
/* based on sim.c: program for running simulations with random and
   deterministic sampling. ML 6/97. */
/*  top-k.c: Michael Littman, Sun Dec 15 19:29:05 EST 1996
 *   Version 4.  Assume entries are positive and sorted (big to small).
 *  Given a set of n integers, list the k subsets that have the
 *  highest sums (in order).
 *
 * Michael Littman, Tue Jun  3 11:38:08 EDT 1997
 *  Modifying to run more standalone.  In particular, does the logit
 *  calculations and sorting itself instead of depending on S to do it.
 * Merlise Clyde, February 2003,  modified to be called from R
 * reworked memory management and tree structures for larger problems
*/

/* Includes. */
#include "sampling.h"

void   update_MCMC_freq(double *MCMC_probs, int *model, int p, int m);
double cond_prob(double *model, int j, int n, double *mean, double *beta_matrix , double eps);

void update_cond_tree(SEXP modelspace, struct Node *tree, SEXP modeldim, struct Var *vars, int p, int n, int kt, int *model, double *real_model, double *marg_probs, double *beta_matrix, double eps);

void  update_Cov(double *Cov, double *priorCov, double *SSgam, double *marg_probs, int n, int m, int print);

void insert_model_tree(struct Node *tree, struct Var *vars,  int n, int *model, int num_models);
 
SEXP amcmc(SEXP Y, SEXP X, SEXP Rprobinit, SEXP Rmodeldim, SEXP incint, SEXP Ralpha,SEXP method, SEXP modelprior, SEXP Rupdate, SEXP Rbestmodel, SEXP Rbestmarg, SEXP plocal, SEXP BURNIN_Iterations, SEXP MCMC_Iterations, SEXP LAMBDA, SEXP DELTA)
{
  SEXP   Rse_m, Rcoef_m, Rmodel_m; 

  SEXP   RXwork = PROTECT(duplicate(X)), RYwork = PROTECT(duplicate(Y));

  int nProtected = 2, nUnique=0, newmodel=0;
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

  double *Xwork, *Ywork, *coefficients,*probs, shrinkage_m, *MCMC_probs,
    SSY, yty, ybar, mse_m, *se_m, MH=0.0, prior_m=1.0, *real_model, prob_i,
    R2_m, RSquareFull, alpha, prone, denom, logmargy, postold, postnew;
  int nobs, p, k, i, j, m, n, l, pmodel, pmodel_old, *xdims, *model_m, *bestmodel, *varin, *varout;
  int mcurrent,  update, n_sure;
  double  mod, rem, problocal, *pigamma, pigammaold, pigammanew, eps, *hyper_parameters;
  double *XtX, *XtY, *XtXwork, *XtYwork, *SSgam, *Cov, *priorCov, *marg_probs;
  double one=1.0, zero=0.0, lambda,  delta, wt = 1.0; 
 
  int inc=1, p2, print = 0;
  int *model, *modelold, bit, *modelwork, old_loc, new_loc;	
  char uplo[] = "U", trans[]="T";
  struct Var *vars;	/* Info about the model variables. */
  NODEPTR tree, branch;

  /* get dimsensions of all variables */


  nobs = LENGTH(Y);
  xdims = INTEGER(getAttrib(X,R_DimSymbol));
  p = xdims[1];
  k = LENGTH(modelprobs);
  update = INTEGER(Rupdate)[0];
  lambda=REAL(LAMBDA)[0];
  delta = REAL(DELTA)[0];
  //  Rprintf("delta %f lambda %f", delta, lambda);
  eps = DBL_EPSILON;
  problocal = REAL(plocal)[0];
  //  Rprintf("Update %i and prob.switch %f\n", update, problocal);
  /* Extract prior on models  */
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));

  /*  Rprintf("n %d p %d \n", nobs, p);  */

  Ywork = REAL(RYwork);
  Xwork = REAL(RXwork);

 
 /* Allocate other variables.  */ 
  XtX  = (double *) R_alloc(p * p, sizeof(double));
  XtXwork  = (double *) R_alloc(p * p, sizeof(double));
  XtY = vecalloc(p);  
  XtYwork = vecalloc(p);


  /* create X matrix */
  for (j=0, l=0; j < p; j++) {
    for (i = 0; i < p; i++) {
      XtX[j*p + i] = 0.0;
    }
    /*    for (i=0; i < nobs; i++) {
       Xmat[i][j] =  REAL(X)[l];
       Xwork[l] = Xmat[i][j];
       l = l + 1; 
       } */
  }
  //  PROTECT(Rprobs = NEW_NUMERIC(p));
  //initprobs = REAL(Rprobinit);


 p2 = p*p;
 ybar = 0.0; SSY = 0.0; yty = 0.0; 

 
 F77_NAME(dsyrk)(uplo, trans, &p, &nobs, &one, &Xwork[0], &nobs, &zero, &XtX[0], &p); 
 yty = F77_NAME(ddot)(&nobs, &Ywork[0], &inc, &Ywork[0], &inc);
 for (i = 0; i< nobs; i++) {
     ybar += Ywork[i];
  }

  ybar = ybar/ (double) nobs;
  SSY = yty - (double) nobs* ybar *ybar;

  F77_NAME(dgemv)(trans, &nobs, &p, &one, &Xwork[0], &nobs, &Ywork[0], &inc, &zero, &XtY[0],&inc);
  
  alpha = REAL(Ralpha)[0];

  vars = (struct Var *) R_alloc(p, sizeof(struct Var));
  probs =  REAL(Rprobs);
  n = sortvars(vars, probs, p); 

  for (i =0; i <n; i++) REAL(MCMCprobs)[vars[i].index] = 0.0;
  MCMC_probs =  REAL(MCMCprobs);


  pigamma = vecalloc(p);
  real_model = vecalloc(n);
  marg_probs = vecalloc(n);
  modelold = ivecalloc(p);
  model = ivecalloc(p);
  modelwork= ivecalloc(p);
  varin= ivecalloc(p);
  varout= ivecalloc(p);


  /* create gamma gamma' matrix */
  SSgam  = (double *) R_alloc(n * n, sizeof(double));
  Cov  = (double *) R_alloc(n * n, sizeof(double));
  priorCov  = (double *) R_alloc(n * n, sizeof(double));
  for (j=0; j < n; j++) {
    for (i = 0; i < n; i++) {
      SSgam[j*n + i] = 0.0;
      Cov[j*n + i] = 0.0;
      priorCov[j*n + i] = 0.0;
      if (j == i)  priorCov[j*n + i] = lambda;
    }
    marg_probs[i] = 0.0;
  }


  /* Make space for the models and working variables. */ 

  /*  pivot = ivecalloc(p); 
  qraux = vecalloc(p);
  work =  vecalloc(2 * p);
  effects = vecalloc(nobs); 
  v =  vecalloc(p * p); 
  betaols = vecalloc(p);
  */

 

  /*  Rprintf("Fit Full Model\n"); */

  if (nobs <= p) {RSquareFull = 1.0;}
  else {
  PROTECT(Rcoef_m = NEW_NUMERIC(p));
  PROTECT(Rse_m = NEW_NUMERIC(p));
  coefficients = REAL(Rcoef_m);  
  se_m = REAL(Rse_m);
  memcpy(coefficients, XtY,  p*sizeof(double));
  memcpy(XtXwork, XtX, p2*sizeof(double));
  memcpy(XtYwork, XtY,  p*sizeof(double));

  mse_m = yty; 
  cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, p, nobs);  

  /*olsreg(Ywork, Xwork,  coefficients, se_m, &mse_m, &p, &nobs, pivot,qraux,work,residuals,effects,v, betaols); */
  RSquareFull =  1.0 - (mse_m * (double) ( nobs - p))/SSY;
  UNPROTECT(2);
  }


  /* fill in the sure things */
  for (i = n, n_sure = 0; i < p; i++)  {
      model[vars[i].index] = (int) vars[i].prob;
      if (model[vars[i].index] == 1) ++n_sure;
  }


  GetRNGstate();
  tree = make_node(-1.0);

  /*  Rprintf("For m=0, Initialize Tree with initial Model\n");  */

  m = 0;
  bestmodel = INTEGER(Rbestmodel);

  INTEGER(modeldim)[m] = n_sure;

  /* Rprintf("Create Tree\n"); */
   branch = tree;

   for (i = 0; i< n; i++) {
      bit =  bestmodel[vars[i].index];
      if (bit == 1) {
	if (i < n-1 && branch->one == NULL) 
	  branch->one = make_node(-1.0);
	if (i == n-1 && branch->one == NULL)
	  branch->one = make_node(0.0);
	branch = branch->one;
      }
      else {
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
  


    /*    Rprintf("Now get model specific calculations \n"); */
 
    pmodel = INTEGER(modeldim)[m];
    PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
    model_m = INTEGER(Rmodel_m);

      for (j = 0, l=0; j < p; j++) {  
	if (model[j] == 1) {
            model_m[l] = j;
           l +=1;}
      }

    SET_ELEMENT(modelspace, m, Rmodel_m);

    Rcoef_m = NEW_NUMERIC(pmodel); PROTECT(Rcoef_m);
    Rse_m = NEW_NUMERIC(pmodel);   PROTECT(Rse_m);
    coefficients = REAL(Rcoef_m);  
    se_m = REAL(Rse_m);

      for (j=0, l=0; j < pmodel; j++) {
        XtYwork[j] = XtY[model_m[j]];
        for  ( i = 0; i < pmodel; i++) {
	  XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
	}
      } 

      
    mse_m = yty; 
    memcpy(coefficients, XtYwork, sizeof(double)*pmodel); 
    cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);  

    R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;

    SET_ELEMENT(beta, m, Rcoef_m);
    SET_ELEMENT(se, m, Rse_m);

    REAL(R2)[m] = R2_m;
    REAL(mse)[m] = mse_m;

    gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
    REAL(sampleprobs)[m] = 1.0;
    REAL(logmarg)[m] = logmargy;
    REAL(shrinkage)[m] = shrinkage_m;
    prior_m = compute_prior_probs(model,pmodel,p, modelprior);
    REAL(priorprobs)[m] = prior_m;
    REAL(Rbestmarg)[0] = REAL(logmarg)[m];
    UNPROTECT(3);


    old_loc = 0;
    pmodel_old = pmodel;
    nUnique=1;
    INTEGER(counts)[0] = 0;
    postold =  REAL(logmarg)[m] + log(REAL(priorprobs)[m]);
    memcpy(modelold, model, sizeof(int)*p);
  /*   Rprintf("model %d max logmarg %lf\n", m, REAL(logmarg)[m]); */

    /*  Rprintf("Now Sample the Rest of the Models \n");  */
    
  
  m = 0;

  //  Need to fix in case the number of sampled models exceeds the space! 
  while (m < INTEGER(BURNIN_Iterations)[0]) {

    memcpy(model, modelold, sizeof(int)*p);
    pmodel =  n_sure;
    MH = 1.0;

    if (pmodel_old == n_sure || pmodel_old == n_sure + n){
	MH =  random_walk(model, vars,  n);
	MH =  1.0 - problocal;
    }
    else {
      if (unif_rand() < problocal) {
      // random
	MH =  random_switch(model, vars, n, pmodel_old, varin, varout );
      }
      else {
      // Randomw walk proposal flip bit//
	MH =  random_walk(model, vars,  n);
      }
    }
    
    branch = tree;
    newmodel= 0;

    for (i = 0; i< n; i++) {
      bit =  model[vars[i].index];
      
      if (bit == 1) {
	if (branch->one != NULL) branch = branch->one;
	else newmodel = 1;
	}
      else {
	if (branch->zero != NULL)  branch = branch->zero;
	else newmodel = 1.0;
      } 
      pmodel  += bit;
    }

    if (pmodel  == n_sure || pmodel == n + n_sure)  MH = 1.0/(1.0 - problocal);

    if (newmodel == 1) {
      new_loc = nUnique;
      PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
      model_m = INTEGER(Rmodel_m);
      for (j = 0, l=0; j < p; j++) {  
	if (model[j] == 1) {
	  model_m[l] = j;
	  l +=1;}
      }	

      Rcoef_m = NEW_NUMERIC(pmodel); PROTECT(Rcoef_m);
      Rse_m = NEW_NUMERIC(pmodel);   PROTECT(Rse_m);
      coefficients = REAL(Rcoef_m);  
      se_m = REAL(Rse_m);
      for (j=0, l=0; j < pmodel; j++) {
        XtYwork[j] = XtY[model_m[j]];
        for  ( i = 0; i < pmodel; i++) {
	  XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
	}	
      }	 

      mse_m = yty; 
      memcpy(coefficients, XtYwork, sizeof(double)*pmodel); 
      cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);  

      R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;
      prior_m = compute_prior_probs(model,pmodel,p, modelprior);
      gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
      postnew = logmargy + log(prior_m);
    }
    else {
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
	SET_ELEMENT(modelspace, nUnique, Rmodel_m);

	SET_ELEMENT(beta, nUnique, Rcoef_m);
	SET_ELEMENT(se, nUnique, Rse_m);

	REAL(R2)[nUnique] = R2_m;
	REAL(mse)[nUnique] = mse_m;
	REAL(sampleprobs)[nUnique] = 1.0;
	REAL(logmarg)[nUnique] = logmargy;
	REAL(shrinkage)[nUnique] = shrinkage_m;
	REAL(priorprobs)[nUnique] = prior_m;
	UNPROTECT(3);
	++nUnique; 
      }

      old_loc = new_loc;
      postold = postnew;
      pmodel_old = pmodel;
      memcpy(modelold, model, sizeof(int)*p);
    }
    else  {
      if (newmodel == 1) UNPROTECT(3);
    }

    INTEGER(counts)[old_loc] += 1;
    
    for (i = 0; i < n; i++) {
     real_model[i] = (double) modelold[vars[i].index];
     REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
   }

   // Update SSgam = gamma gamma^T + SSgam 
   F77_NAME(dsyr)("U", &n,  &one, &real_model[0], &inc,  &SSgam[0], &n);
   m++;
  }
  
  //  Rprintf("\n%d Unique models sampled during burnin\n", nUnique);


// Compute marginal probabilities  
  mcurrent = nUnique;
  compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
  compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);        
  for (i = 0; i < n; i++) {
    marg_probs[i] = wt*(REAL(MCMCprobs)[vars[i].index]/ (double) m) + 
      (1.0 - wt)* probs[vars[i].index];
  }	
  //  print=1;
  update_Cov(Cov, priorCov, SSgam, marg_probs, n, m, print);
 
// Global-Proposal
// Initialize post old proposal
 pigammaold = 0.0;
 for (i = 0; i < n; i++) {
   if (modelold[vars[i].index] == 1 ){	
     real_model[i] = 1.0;
     pigammaold += log(cond_prob(real_model,i, n, marg_probs,Cov, delta));
   }
   else {
     real_model[i] = 0.0;
     pigammaold += log(1.0 - cond_prob(real_model,i, n, marg_probs,Cov, delta));
   }
 }
 
//  need to fix to make sure that nUnique is less than nModels 

 while (m < INTEGER(BURNIN_Iterations)[0] + INTEGER(MCMC_Iterations)[0]) {  
      // for (m = 0; m < k; m++) {
   memcpy(model, modelold, sizeof(int)*p);
   pmodel =  n_sure;
   MH = 1.0;
   pigammanew = 0.0;
   branch = tree;
   newmodel = 0;
   for (i = 0; i < n; i++) {
     prob_i = cond_prob(real_model,i, n, marg_probs,Cov,delta);
     bit =  withprob(prob_i);
     
     if (bit == 1) {
       pigammanew += log(prob_i);
       if (branch->one != NULL) branch = branch->one;
       else newmodel= 1;
     }
     else {
       pigammanew += log(1.0 - prob_i);
       if (branch->zero != NULL) branch = branch->zero;
       else newmodel= 1;
     } 
     model[vars[i].index] = bit; 
     real_model[i] = (double) bit;
     pmodel  += bit;
   }
   if (newmodel == 1) {
     new_loc = nUnique;
     PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
     model_m = INTEGER(Rmodel_m);
     for (j = 0, l=0; j < p; j++) {  
       if (model[j] == 1) {
	 model_m[l] = j;
	 l +=1;}
     }	

     Rcoef_m = NEW_NUMERIC(pmodel); PROTECT(Rcoef_m);
     Rse_m = NEW_NUMERIC(pmodel);   PROTECT(Rse_m);
     coefficients = REAL(Rcoef_m);  
     se_m = REAL(Rse_m);
     for (j=0, l=0; j < pmodel; j++) {
       XtYwork[j] = XtY[model_m[j]];
       for  ( i = 0; i < pmodel; i++) {
	 XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
       }		
      }	 

     mse_m = yty; 
     memcpy(coefficients, XtYwork, sizeof(double)*pmodel); 
     cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);  
     
     R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;
     prior_m = compute_prior_probs(model,pmodel,p, modelprior);
     gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
     postnew = logmargy + log(prior_m);
    }	
   else {
     new_loc = branch->where;
     postnew =  REAL(logmarg)[new_loc] + log(REAL(priorprobs)[new_loc]);      
   } 

   MH = exp(postnew - postold + pigammaold - pigammanew);
   // Rprintf("it %d MH %lf  new %lf old %lf propold %lf propnew %lf \n", m, MH, postnew, postold, pigammanew, pigammaold);
   if (unif_rand() < MH) {
     
     if (newmodel ==1)  {
       new_loc = nUnique;

       insert_model_tree(tree, vars, n, model, nUnique);

       INTEGER(modeldim)[nUnique] = pmodel;
       SET_ELEMENT(modelspace, nUnique, Rmodel_m);
       
       SET_ELEMENT(beta, nUnique, Rcoef_m);
       SET_ELEMENT(se, nUnique, Rse_m);
       
       REAL(R2)[nUnique] = R2_m;
       REAL(mse)[nUnique] = mse_m;
       
       REAL(logmarg)[nUnique] = logmargy;
       REAL(shrinkage)[nUnique] = shrinkage_m;
       REAL(priorprobs)[nUnique] = prior_m;
       UNPROTECT(3);
       ++nUnique; 
     }	

     old_loc = new_loc;
     pigammaold = pigammanew;
     REAL(sampleprobs)[old_loc] = pigammaold;
     postold = postnew;
     pmodel_old = pmodel;
     memcpy(modelold, model, sizeof(int)*p);
   }
   else  {
     if (newmodel == 1) UNPROTECT(3);
   }	
 
   INTEGER(counts)[old_loc] += 1;

   for (i = 0; i < n; i++) {
     real_model[i] = (double) modelold[vars[i].index];
     REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
   }
   F77_NAME(dsyr)("U", &n,  &one, &real_model[0], &inc,  &SSgam[0], &n);
   m++;

   
   rem = modf((double) m/(double) update, &mod);
   if (rem  == 0.0) {
     mcurrent = nUnique;
     compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
     compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);        

     for (i = 0; i < n; i++) {
       marg_probs[i] = wt*(REAL(MCMCprobs)[vars[i].index]/ (double) m) + 
	 (1.0 - wt)*probs[vars[i].index];
     }
     update_Cov(Cov, priorCov, SSgam, marg_probs, n, m, print);
     // Initialize post old proposal
     pigammaold = 0.0;
     for (i = 0; i < n; i++) {
       if (modelold[vars[i].index] == 1 ){	
	 real_model[i] = 1.0;
	 pigammaold += log(cond_prob(real_model,i, n, marg_probs,Cov, delta));
       }
       else {
	 real_model[i] = 0.0;
	 pigammaold += log(1.0 - cond_prob(real_model,i, n, marg_probs,Cov, delta));
       }}
   }
 }

 mcurrent = nUnique;
 compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
 compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);        
 wt = 0.1;
 for (i = 0; i < n; i++) {
   marg_probs[i] = wt* (REAL(MCMCprobs)[vars[i].index]/ (double) m) + 
				    (1.0 - wt)*probs[vars[i].index];
   //  marg_probs[n-1-i] =  REAL(MCMCprobs)[vars[i].index]/ (double) m;
  REAL(MCMCprobs)[vars[i].index] /= (double) m;
  }

 update_Cov(Cov, priorCov, SSgam, marg_probs, n, m, print);

//  Now sample W/O Replacement 
// Rprintf("NumUnique Models Accepted %d \n", nUnique);
 INTEGER(NumUnique)[0] = nUnique;

 if (nUnique < k) {

    //    compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
    //    compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);        
    
   //   update_MCMC_probs(MCMC_probs, vars, n, p);
   //   Rprintf("Update Tree\n");

   update_cond_tree(modelspace, tree, modeldim, vars, p, n, nUnique, modelwork, real_model, marg_probs, Cov, eps);      
  
         
   //   Rprintf("\nNow sample the rest without replacement\n");
   
   for (m = nUnique; m < k; m++) {
     INTEGER(modeldim)[m]  = n_sure;
     
     branch = tree;
  
     for (i = 0; i< n; i++) {
       pigamma[i] = 1.0;

       if (branch->prob == -1.0) {
	 branch->prob = cond_prob(real_model,i, n, marg_probs,Cov, delta);
	 branch->update = m+mcurrent;
       }
       	 
       bit =  withprob(branch->prob);
       real_model[n-i-1] = (double) bit;
     
       if (bit == 1) {
	 for (j=0; j<=i; j++)  pigamma[j] *= branch->prob;
	 if (i < n-1 && branch->one == NULL) 
	 //	    branch->one = make_node(vars[i+1].prob);
	 branch->one = make_node(cond_prob(real_model,i+1, n, marg_probs,Cov , delta));
       if (i == n-1 && branch->one == NULL)
	 branch->one = make_node(0.0);
       branch = branch->one;
       }
     else {
       for (j=0; j<=i; j++)  pigamma[j] *= (1.0 - branch->prob);
       if (i < n-1 && branch->zero == NULL)
	 //	 branch->zero = make_node(vars[i+1].prob);
	 branch->zero = make_node(cond_prob(real_model,i+1, n, marg_probs,Cov, delta));
       if (i == n-1 && branch->zero == NULL)
	 branch->zero = make_node(0.0);
       branch = branch->zero;
     }
       model[vars[i].index] = bit; 
       INTEGER(modeldim)[m]  += bit;
     }

     REAL(sampleprobs)[m] = pigamma[0]; 

     pmodel = INTEGER(modeldim)[m];

    // Now subtract off the visited probability mass. 
     branch=tree;
     for (i = 0; i < n; i++) {
       bit = model[vars[i].index];
       prone = branch->prob;
       if (bit == 1) prone -= pigamma[i];
       denom = 1.0 - pigamma[i];
       if (denom <= 0.0) {
	 if (denom < 0.0) {
	   Rprintf("neg denominator %le %le %le !!!\n", pigamma, denom, prone);
	   if (branch->prob < 0.0 && branch->prob < 1.0)
	     Rprintf("non extreme %le\n", branch->prob);}
	 denom = 0.0;}
       else {
	 if  (prone <= 0)  prone = 0.0;
	 if  (prone > denom)  {
	   if (prone <= eps) prone = 0.0;
	   else prone = 1.0;
	  // Rprintf("prone > 1 %le %le %le %le !!!\n", pigamma, denom, prone, eps);
	 }
	 else prone = prone/denom;
       }
       if (prone > 1.0 || prone < 0.0) 
	 	 Rprintf("%d %d Probability > 1!!! %le %le  %le %le \n",
	 	 m, i, prone, branch->prob, denom, pigamma);
      //      if (bit == 1)  pigamma /= (branch->prob);
      //	      else  pigamma /= (1.0 - branch->prob); 
      //	      if (pigamma > 1.0) pigamma = 1.0; 
       branch->prob  = prone;
       if (bit == 1) branch = branch->one;
       else  branch = branch->zero;
      //      Rprintf("%d %d \n",  branch->done, n - i); 
      //      if (log((double) branch->done) < (n - i)*log(2.0)) {
      //	if (bit == 1) branch = branch->one;
      //	else  branch = branch->zero;
      //}
      //else {
      //	    branch->one = NULL;
      //	    branch->zero = NULL; 
      //	    break; } 
     }
    /* Now get model specific calculations */ 

    PROTECT(Rmodel_m = allocVector(INTSXP, pmodel));
    model_m = INTEGER(Rmodel_m);

    for (j = 0, l=0; j < p; j++) {  
      if (model[j] == 1) {
	model_m[l] = j;
	l +=1;}
      }	
 

      SET_ELEMENT(modelspace, m, Rmodel_m);
   
      for (j=0, l=0; j < pmodel; j++) {
	XtYwork[j] = XtY[model_m[j]];
        for  ( i = 0; i < pmodel; i++) {
	  XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
	}	
      } 

    
      PROTECT(Rcoef_m = allocVector(REALSXP,pmodel));
      PROTECT(Rse_m = allocVector(REALSXP,pmodel));
      coefficients = REAL(Rcoef_m);  
      se_m = REAL(Rse_m);
  
      mse_m = yty; 
      memcpy(coefficients, XtYwork, sizeof(double)*pmodel); 
      cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);  

  
    //    olsreg(Ywork, Xwork, coefficients, se_m, &mse_m, &pmodel, &nobs, pivot,qraux,work,residuals,effects,v,betaols);   
  
      R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;

      SET_ELEMENT(beta, m, Rcoef_m);
      SET_ELEMENT(se, m, Rse_m);

      REAL(R2)[m] = R2_m;
      REAL(mse)[m] = mse_m;

   gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0],  RSquareFull, SSY, &logmargy, &shrinkage_m);
   REAL(logmarg)[m] = logmargy;
   REAL(shrinkage)[m] = shrinkage_m;
   REAL(priorprobs)[m] = compute_prior_probs(model,pmodel,p, modelprior);
   //   REAL(priorprobs)[m] = beta_binomial(pmodel,p, hyper_parameters);
   //  if (REAL(logmarg)[m] > REAL(Rbestmarg)[0]) {
   //   for (i=0; i < p; i++) {
   //	 bestmodel[i] = model[i];}
   //       REAL(Rbestmarg)[0] = REAL(logmarg)[m];
   //} 
       
   UNPROTECT(3);  
   }	
 }


 // Rprintf("modelprobs\n");
 compute_modelprobs(modelprobs, logmarg, priorprobs,k);
 // Rprintf("marginal probs\n");
 compute_margprobs(modelspace, modeldim, modelprobs, probs, k, p);  
 
 // Rprintf("saving\n");
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

  SET_VECTOR_ELT(ANS, 6, mse);
  SET_STRING_ELT(ANS_names, 6, mkChar("mse"));

  SET_VECTOR_ELT(ANS, 7, beta);
  SET_STRING_ELT(ANS_names, 7, mkChar("mle"));

  SET_VECTOR_ELT(ANS, 8, se);
  SET_STRING_ELT(ANS_names, 8, mkChar("mle.se"));

  SET_VECTOR_ELT(ANS, 9, shrinkage);
  SET_STRING_ELT(ANS_names, 9, mkChar("shrinkage"));

  SET_VECTOR_ELT(ANS, 10, modeldim);
  SET_STRING_ELT(ANS_names, 10, mkChar("size"));
 
  SET_VECTOR_ELT(ANS, 11, R2);
  SET_STRING_ELT(ANS_names, 11, mkChar("R2"));

  SET_VECTOR_ELT(ANS, 12, counts);
  SET_STRING_ELT(ANS_names, 12, mkChar("freq"));

  SET_VECTOR_ELT(ANS, 13, MCMCprobs);
  SET_STRING_ELT(ANS_names, 13, mkChar("probs.MCMC"));

  SET_VECTOR_ELT(ANS, 14, NumUnique);
  SET_STRING_ELT(ANS_names, 14, mkChar("n.Unique"));

  setAttrib(ANS, R_NamesSymbol, ANS_names);
  UNPROTECT(nProtected);
  //  Rprintf("Return\n");
  PutRNGstate();

  return(ANS);  
}

double random_walk(int *model, struct Var *vars, int n) {
  int index;
      index = ftrunc(n*unif_rand());
      model[vars[index].index] = 1 - model[vars[index].index];
      return(1.0);
 }

double random_switch(int *model, struct Var *vars, int n, int pmodel, int *varin, int *varout) {
  int  j, k, swapin, swapout, num_to_swap_in, num_to_swap_out;


      j = 0; k = 0;
      while (j < n && k < pmodel)
        {
	  if (model[vars[j].index]==1) {varin[k] = vars[j].index; k++;}
	  j++ ;
	}
      num_to_swap_in = k;
      j = 0; k = 0;

      while (j< n)
        {
	  if (model[vars[j].index]==0) {varout[k] = vars[j].index; k++;}
	  j++ ;
        }	
      num_to_swap_out = k;	

      swapin = ftrunc(unif_rand()*num_to_swap_in);    // swapin :corresponds to position of randomly chosen included variable
      swapout = ftrunc(unif_rand()*num_to_swap_out);  // swapout :corresponds to position of randomly chosen excluded variable

      model[varin[swapin]] = 0;
      model[varout[swapout]] =1; 
   

    return(1.0);
}


double cond_prob(double *model, int j, int n, double *mean, double *beta_matrix , double delta) {
  double  prob;
  int i;
 
  prob = mean[j]; 
  //  Rprintf("%f\n", prob);
  for (i = 0; i < j; i++) {
    prob += -beta_matrix[j*n + i]*(model[i] - mean[i])/beta_matrix[j*n+j];
    //    Rprintf("model %f beta %f mean %f ", model[i], beta_matrix[i*n + j], mean[i]);
  }
  //  Rprintf("\n%f ", prob);
  if (prob <= 0.0) prob = .02;
  if (prob >= 1.0) prob = 1.0 - .02;
  //    Rprintf("%f \n", prob);
  return(prob);
}

void update_cond_tree(SEXP modelspace, struct Node *tree, SEXP modeldim, struct Var *vars, int p, int n, int kt, int *model, double *real_model, double *marg_probs, double *beta_matrix, double delta)
{
  int i,m, bit;
  double prone, pigamma, przero;
  SEXP model_m;
  struct Node *branch;


  for (m = 0;  m <= kt; m++) {
    branch = tree;
    PROTECT(model_m = VECTOR_ELT(modelspace, m));

    for (i = 0; i < p; i++) model[i] = 0;
    for (i = 0; i < INTEGER(modeldim)[m]; i++) 
       model[INTEGER(model_m)[i]] = 1;
    for (i = 0; i < p; i++) real_model[p-1-i] = (double) model[vars[i].index];

    pigamma = 0.0;
    for (i = 0; i < n; i++) {
  	if (branch->update != kt) {
	  //          branch->prob = vars[i].prob;
	  branch->prob = cond_prob(real_model,i,n, marg_probs, beta_matrix, delta);
          branch->update = kt;
	}
  	bit = model[vars[i].index];
        if (bit ==  1)  {
	  pigamma += log(branch->prob);
	  branch = branch->one;}
        else{
	  pigamma += log(1.0 - branch->prob);
	  branch = branch->zero;}
    }

    branch = tree;
    for (i = 0; i < n; i++) {
      bit = model[vars[i].index];
      if (bit == 1) {
	prone = (branch->prob - exp(pigamma));
        przero = 1.0 - branch->prob;
        pigamma -= log(branch->prob);}
      else {
	prone = branch->prob;
        przero = 1.0 - branch->prob  - exp(pigamma);
        pigamma -= log(1.0 - branch->prob);}
      if  (prone <= 0.0 )  prone = 0.;
      if  (przero <= 0.0 )  przero = 0.;
      branch->prob  = prone/(prone + przero);
      if (prone <= 0.0) branch->prob = 0.;

      if (bit == 1) branch = branch->one;
      else branch = branch->zero;
    }
    UNPROTECT(1);
  }
}



void  update_Cov(double *Cov, double *priorCov, double *SSgam, double *marg_probs, int n, int m, int print) {
  double alpha, one=1.0;
  int matsize=n*n, inc=1, i,j,info=1;
    
  memcpy(Cov, SSgam, sizeof(double)*matsize);
  if (print == 1) {
   Rprintf("SS: %d iterations\n", m);

    for (j=0; j < n; j++ ){
      Rprintf("%d  %f  ", j, marg_probs[j]);
      for(i = 0; i < n; i++){
	Rprintf("%f ", Cov[i*n + j]);}
      Rprintf("\n");
    }
  }
  alpha  = -(double) m;
  // Compute SSgam - m prob prob^T 
  F77_NAME(dsyr)("U", &n,  &alpha, &marg_probs[0], &inc,  &Cov[0], &n);
  alpha = 1.0/(double) m;
  // divide SSgam - m prob prob^T by m -> Cov
  F77_NAME(dscal)(&matsize,  &alpha,  &Cov[0], &inc);
  // Add priorCov to Cov
  F77_NAME(daxpy)(&matsize, &one, &priorCov[0], &inc, &Cov[0], &inc);
  
  if (print == 1) {
    Rprintf("Cov:\n");

    for (j=0; j < n; j++ ){
      for(i = 0; i < n; i++){
	Rprintf("%f ", Cov[i*n + j]);}
      Rprintf("\n");
    }
  }
  F77_NAME(dpotrf)("U", &n,  &Cov[0], &n, &info);
  F77_NAME(dtrtri)("U","N", &n, &Cov[0], &n, &info );

  if (print == 1) {
    Rprintf("inverse of Chol(Cov(SSgam)):\n");
    for (j=0; j < n; j++ ){
      for(i = 0; i < n; i++){
	Rprintf("%f ", Cov[i*n + j]);}
      Rprintf("\n");
    } 
  }
}
