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
#include "bas.h"

void   update_MCMC_freq(double *MCMC_probs, int *model, int p, int m);
double cond_prob(double *model, int j, int n, double *mean, double *beta_matrix , double eps);

void update_cond_tree(SEXP modelspace, struct Node *tree, SEXP modeldim, struct Var *vars, int p, int n, int kt, int *model, double *real_model, double *marg_probs, double *beta_matrix, double eps);

void  update_Cov(double *Cov, double *priorCov, double *SSgam, double *marg_probs, int n, int m, int print);

void insert_model_tree(struct Node *tree, struct Var *vars,  int n, int *model, int num_models);

// [[register]]
SEXP mcmcbas(SEXP Y, SEXP X, SEXP Rweights, SEXP Rprobinit, SEXP Rmodeldim, 
             SEXP incint, SEXP Ralpha,SEXP method, SEXP modelprior, SEXP Rupdate, 
             SEXP Rbestmodel,  SEXP plocal,
             SEXP BURNIN_Iterations, SEXP MCMC_Iterations, SEXP LAMBDA, SEXP DELTA,
             SEXP Rparents)
{

  int nProtected = 0;
  SEXP RXwork = PROTECT(duplicate(X)); nProtected++; 
  SEXP RYwork = PROTECT(duplicate(Y));  nProtected++; 
  int nModels=LENGTH(Rmodeldim);
  int nUnique=0, newmodel=0;
  
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
  SEXP Rse_m = NULL, Rcoef_m = NULL, Rmodel_m;

  double *Xwork, *Ywork, *wts, *coefficients,*probs, shrinkage_m,
         mse_m, *se_m, MH=0.0, prior_m=1.0, *real_model,
         R2_m, RSquareFull, prone, denom, logmargy, postold, postnew;
  int  i, j, m, n, l, pmodel_old, *model_m, *bestmodel, *varin, *varout;
  int mcurrent,  update, n_sure;
  double  mod, rem, problocal,  eps;
//  double *hyper_parameters, delta;  //For future use
  // double *marg_probs;
 // double *SSgam, *Cov, *priorCov;

  
  int *modelold, bit, *modelwork, old_loc, new_loc;
  //  char uplo[] = "U", trans[]="T";
  struct Var *vars;	/* Info about the model variables. */
  NODEPTR tree, branch;

  /* get dimsensions of all variables */
  int nobs = LENGTH(Y);
  int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
  int k = LENGTH(modelprobs);
 
  // double lambda=REAL(LAMBDA)[0];
  double alpha = REAL(Ralpha)[0];
  SEXP Rbestmodel_new = PROTECT(duplicate(Rbestmodel)); nProtected++;
  
  
  update = INTEGER(Rupdate)[0];
  eps = DBL_EPSILON;
  problocal = REAL(plocal)[0];
  
 // hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));


  Ywork = REAL(RYwork);
  Xwork = REAL(RXwork);
  wts = REAL(Rweights);


 /* Allocate other variables.  */
  double *XtX, *XtY, *XtXwork, *XtYwork, yty, SSY;
  PrecomputeData(Xwork, Ywork, wts, &XtXwork, &XtYwork, &XtX, &XtY, &yty, &SSY, p, nobs);

  

  vars = (struct Var *) R_alloc(p, sizeof(struct Var));
  probs =  REAL(Rprobs);
  n = sortvars(vars, probs, p);

  for (i =n; i <p; i++) REAL(MCMCprobs)[vars[i].index] = probs[vars[i].index];
  for (i =0; i <n; i++) REAL(MCMCprobs)[vars[i].index] = 0.0;

  int noInclusionIs1 = no_prior_inclusion_is_1(p, probs);

  

  real_model = vecalloc(n);
//  marg_probs = vecalloc(n);
  modelold = ivecalloc(p);
  modelwork= ivecalloc(p);
  varin= ivecalloc(p);
  varout= ivecalloc(p);


  
  //  allocate working model and fill in the sure things
  int *model = ivecalloc(p);
  memset(model, 0, p * sizeof(int));
  
  

  /* fill in the sure things */
  for (i = n, n_sure = 0; i < p; i++)  {
      model[vars[i].index] = (int) vars[i].prob;
      if (model[vars[i].index] == 1) ++n_sure;
  }

  RSquareFull = CalculateRSquareFull(XtY, XtX, XtXwork, XtYwork, Rcoef_m, Rse_m, p, nobs, yty, SSY);
   
  GetRNGstate();
  tree = make_node(-1.0);

  /*  Rprintf("For m=0, Initialize Tree with initial Model\n");  */

  m = 0;
  bestmodel = INTEGER(Rbestmodel_new);
  INTEGER(modeldim)[m] = n_sure;

  /* Rprintf("Create Tree\n"); */

   double *pigamma = vecalloc(p);
   branch = tree;
   CreateTree_with_pigamma(branch, vars, bestmodel, model, n, m,
                           modeldim, pigamma, Rparents);
   
   branch=tree;
   Substract_visited_probability_mass(branch, vars, model, n, m, pigamma,eps);
   

  //  Rprintf("Now get model specific calculations \n"); 

    int pmodel = INTEGER(modeldim)[m];
    PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
    PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
    PROTECT(Rse_m = NEW_NUMERIC(pmodel));
    
    model_m = GetModel_m(Rmodel_m, model, p);

    SET_ELEMENT(modelspace, m, Rmodel_m);

    coefficients = REAL(Rcoef_m);
    se_m = REAL(Rse_m);

      for (j=0, l=0; j < pmodel; j++) {
        XtYwork[j] = XtY[model_m[j]];
        for  ( i = 0; i < pmodel; i++) {
	         XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
	}
      }

    R2_m = 0.0;
    mse_m = yty;
    memcpy(coefficients, XtYwork, sizeof(double)*pmodel);
    cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);

    if (pmodel > 1)   R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;
    
    
    SET_ELEMENT(beta, m, Rcoef_m);
    SET_ELEMENT(se, m, Rse_m);

    REAL(R2)[m] = R2_m;
    REAL(mse)[m] = mse_m;

    gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);

    REAL(sampleprobs)[m] = 1.0;
    REAL(logmarg)[m] = logmargy;
    REAL(shrinkage)[m] = shrinkage_m;
    prior_m  = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
    REAL(priorprobs)[m] = prior_m;

    UNPROTECT(3);

    
    old_loc = 0;
    pmodel_old = pmodel;
    nUnique=1;
    INTEGER(counts)[0] = 1;
    postold =  REAL(logmarg)[m] + log(REAL(priorprobs)[m]);
    memcpy(modelold, model, sizeof(int)*p);
  /*   Rprintf("model %d max logmarg %lf\n", m, REAL(logmarg)[m]); */

  Rprintf("Now Sample Models with MCMC  \n");  


  m = 0;

  while (nUnique < k && m < INTEGER(BURNIN_Iterations)[0]) {
  // Rprintf("m = %d \n", m);
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
      R2_m = 0.0;
      mse_m = yty;
      memcpy(coefficients, XtYwork, sizeof(double)*pmodel);
      cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);
      if (pmodel > 1)  R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;
      prior_m = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
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
      /* store in opposite order so nth variable is first */
     real_model[n-1-i] = (double) modelold[vars[i].index];
     REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
   }

   m++;
  }

 for (i = 0; i < n; i++) {
     REAL(MCMCprobs)[vars[i].index] /= (double) m;
 }
  Rprintf("\nNum Unique%d  in %d MCMC iterations \n", nUnique, m);


// Compute marginal probabilities
  mcurrent = nUnique;
  compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
  compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);



//  Now sample W/O Replacement
 Rprintf("Sample w/out Replacement Now \n", nUnique);
 
 INTEGER(NumUnique)[0] = nUnique;


 if (nUnique < k) {
   update_probs(probs, vars, mcurrent, k, p);
   Rprintf("updating tree for SWOR\n");
   update_tree(modelspace, tree, modeldim, vars, k,p,n,mcurrent, modelwork);
   Rprintf("Done!\n");
  for (m = nUnique;  m < k; m++) {
    for (i = n; i < p; i++)  {
      INTEGER(modeldim)[m]  +=  model[vars[i].index];
    }

    branch = tree;

    for (i = 0; i< n; i++) {
      pigamma[i] = 1.0;
      bit =  withprob(branch->prob);

      /*    branch->done += 1; */

	if (bit == 1) {
	  for (j=0; j<=i; j++)  pigamma[j] *= branch->prob;
	  if (i < n-1 && branch->one == NULL)
	    branch->one = make_node(vars[i+1].prob);
          if (i == n-1 && branch->one == NULL)
	    branch->one = make_node(0.0);
	  branch = branch->one;
	}
        else {
	  for (j=0; j<=i; j++)  pigamma[j] *= (1.0 - branch->prob);
	  if (i < n-1 && branch->zero == NULL)
	    branch->zero = make_node(vars[i+1].prob);
          if (i == n-1 && branch->zero == NULL)
	    branch->zero = make_node(0.0);
	  branch = branch->zero;
	  }
	model[vars[i].index] = bit;
	INTEGER(modeldim)[m]  += bit;
    }

    REAL(sampleprobs)[m] = pigamma[0];
    pmodel = INTEGER(modeldim)[m];

    /* Now subtract off the visited probability mass. */
    branch=tree;
    for (i = 0; i < n; i++) {
      bit = model[vars[i].index];
      prone = branch->prob;
      if (bit == 1) prone -= pigamma[i];
      denom = 1.0 - pigamma[i];
      if (denom <= 0.0) {
	if (denom < 0.0) {
	  warning("neg denominator %le %le %le !!!\n", pigamma, denom, prone);
	  if (branch->prob < 0.0 && branch->prob < 1.0)
	    warning("non extreme %le\n", branch->prob);}
        denom = 0.0;}
      else {
	if  (prone <= 0)  prone = 0.0;
	if  (prone > denom)  {
          if (prone <= eps) prone = 0.0;
	  else prone = 1.0;
	  /* Rprintf("prone > 1 %le %le %le %le !!!\n", pigamma, denom, prone, eps);*/
	}
	else prone = prone/denom;
      }
      if (prone > 1.0 || prone < 0.0)
		Rprintf("%d %d Probability > 1!!! %le %le  %le %le \n",
		m, i, prone, branch->prob, denom, pigamma);


      /*      if (bit == 1)  pigamma /= (branch->prob);
	      else  pigamma /= (1.0 - branch->prob);
	      if (pigamma > 1.0) pigamma = 1.0; */
      branch->prob  = prone;
      if (bit == 1) branch = branch->one;
      else  branch = branch->zero;
    }

    /* Now get model specific calculations */

    PROTECT(Rmodel_m = allocVector(INTSXP, pmodel));
    PROTECT(Rcoef_m = allocVector(REALSXP,pmodel));
    PROTECT(Rse_m = allocVector(REALSXP,pmodel));
    coefficients = REAL(Rcoef_m);
    se_m = REAL(Rse_m);
    model_m = INTEGER(Rmodel_m);

  for (j = 0, l=0; j < p; j++) {
	    if (model[j] == 1) {
          model_m[l] = j;
           l +=1;
           }
  }

     SET_ELEMENT(modelspace, m, Rmodel_m);

      for (j=0, l=0; j < pmodel; j++) {
         XtYwork[j] = XtY[model_m[j]];
        for  ( i = 0; i < pmodel; i++) {
        	 XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
	      }
      }

    mse_m = yty;
    memcpy(coefficients, XtYwork, sizeof(double)*pmodel);
    cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);


/*    olsreg(Ywork, Xwork, coefficients, se_m, &mse_m, &pmodel, &nobs, pivot,qraux,work,residuals,effects,v,betaols);   */
    if (pmodel > 1)  R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;

    SET_ELEMENT(beta, m, Rcoef_m);
    SET_ELEMENT(se, m, Rse_m);

    REAL(R2)[m] = R2_m;
    REAL(mse)[m] = mse_m;

   gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0],  RSquareFull, SSY, &logmargy, &shrinkage_m);
   REAL(logmarg)[m] = logmargy;
   REAL(shrinkage)[m] = shrinkage_m;
   REAL(priorprobs)[m] = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);


    if (m > 1) {
      rem = modf((double) m/(double) update, &mod);
      if (rem  == 0.0) {
      	mcurrent = m;
	      compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
	      compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);
	      if (update_probs(probs, vars, mcurrent, k, p) == 1) {
        	  Rprintf("Updating Model Tree %d \n", m);
	           update_tree(modelspace, tree, modeldim, vars, k,p,n,mcurrent, modelwork);
	      }
      }}
  UNPROTECT(3);
  }
 }

 compute_modelprobs(modelprobs, logmarg, priorprobs,k);
 compute_margprobs(modelspace, modeldim, modelprobs, probs, k, p);

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




