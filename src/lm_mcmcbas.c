
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
             SEXP Rthin, SEXP Rparents, SEXP Rpivot, SEXP Rtol)
{

  int nProtected = 0;
  SEXP RXwork = PROTECT(duplicate(X)); nProtected++; 
  SEXP RYwork = PROTECT(duplicate(Y));  nProtected++; 
  int nModels=LENGTH(Rmodeldim);
  int pivot = LOGICAL(Rpivot)[0];
  double tol = REAL(Rtol)[0];
  int nUnique=0;
  
  //  Rprintf("Allocating Space for %d Models\n", nModels) ;
  SEXP ANS = PROTECT(allocVector(VECSXP, 16)); ++nProtected;
  SEXP ANS_names = PROTECT(allocVector(STRSXP, 16)); ++nProtected;
  SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
  SEXP MCMCprobs= PROTECT(duplicate(Rprobinit)); ++nProtected;
  SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
  SEXP rank = PROTECT(allocVector(INTSXP, nModels)); ++nProtected;
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

  double *Xwork, *Ywork, *wts, *probs, shrinkage_m,
         mse_m,  MH=0.0, 
         R2_m, RSquareFull, logmargy, postold, postnew;
  int  i, m, n, pmodel_old, *model_m, *bestmodel, rank_m;
  int mcurrent, n_sure;
  double  mod, rem;
    //  char uplo[] = "U", trans[]="T";
 

  /* get dimsensions of all variables */
  int nobs = LENGTH(Y);
  int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
  int k = LENGTH(modelprobs);
  double alpha = REAL(Ralpha)[0];
  int thin = INTEGER(Rthin)[0];
  SEXP Rbestmodel_new = PROTECT(duplicate(Rbestmodel)); nProtected++;
  
  
  int update = INTEGER(Rupdate)[0];
  double eps = DBL_EPSILON;
  
  
 // hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));


  Ywork = REAL(RYwork);
  Xwork = REAL(RXwork);
  wts = REAL(Rweights);


  double *XtX, *XtY, *XtXwork, *XtYwork, yty= 0.0, SSY=0.0;
  PrecomputeData(Xwork, Ywork, wts, &XtXwork, &XtYwork, &XtX, &XtY, &yty, &SSY, p, nobs);

  struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var));
  probs =  REAL(Rprobs);
  n = sortvars(vars, probs, p);
  int noInclusionIs1 = no_prior_inclusion_is_1(p, probs);

  for (i =n; i <p; i++) REAL(MCMCprobs)[vars[i].index] = probs[vars[i].index];
  for (i =0; i <n; i++) REAL(MCMCprobs)[vars[i].index] = 0.0;

  
  
  SEXP Rse_m = NULL, Rcoef_m = NULL, Rmodel_m=NULL;
  RSquareFull = CalculateRSquareFull(XtY, XtX, XtXwork, XtYwork, Rcoef_m, Rse_m, p, nobs, yty, SSY);
   
   //  allocate working model and fill in the sure things
   int *model = ivecalloc(p);
   memset(model, 0, p * sizeof(int));
   
   int   *modelwork = ivecalloc(p);
   memset(modelwork, 0, p * sizeof(int));
   
   
   /* fill in the sure things */
   for (i = n, n_sure = 0; i < p; i++)  {
     model[vars[i].index] = (int) vars[i].prob;
     if (model[vars[i].index] == 1) ++n_sure;
   }
   
  
  GetRNGstate();
  
  NODEPTR tree, branch;
  tree = make_node(vars[0].prob);
  

  /*  Rprintf("For m=0, Initialize Tree with initial Model\n");  */

  m = 0;
  bestmodel = INTEGER(Rbestmodel_new);
  REAL(logmarg)[m] = 0.0;
  INTEGER(modeldim)[m] = 0;
  
  for (i = n; i < p; i++)  {
    model[vars[i].index] = bestmodel[vars[i].index];
    INTEGER(modeldim)[m]  +=  bestmodel[vars[i].index];
  }

  /* Rprintf("Create Tree\n"); */
  
   double *pigamma = vecalloc(p);
   branch = tree;
   CreateTree_with_pigamma(branch, vars, bestmodel, model, n, m,
                           modeldim, pigamma, Rparents);
   
   Rprintf("initialized tree\n");
   
   branch=tree;
   Substract_visited_probability_mass(branch, vars, model, n, m, pigamma,eps);
   
   Rprintf("updated prob tree\n");
   
  //  Rprintf("Now get model specific calculations \n"); 

    int pmodel = INTEGER(modeldim)[m];
    PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
    memset(INTEGER(Rmodel_m), 0, pmodel * sizeof(int));
    PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
    PROTECT(Rse_m = NEW_NUMERIC(pmodel));
    
    model_m = GetModel_m(Rmodel_m, model, p);
    
    R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel,
                    p, nobs, m, &mse_m, &rank_m, pivot, tol);
    INTEGER(rank)[0] = rank_m;
    
    gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy,
                  &shrinkage_m);
    
    
    double prior_m  = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
    if (prior_m == 0.0)  warning("warning initial model has 0 prior probabilty\n");
  
    SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
    SetModel(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2, m);
    
    
    int newmodel=0, nsamples=0;
    double *real_model = vecalloc(n);
    int *modelold = ivecalloc(p);
    int old_loc = 0;
    int new_loc;
    pmodel_old = pmodel;
    nUnique=1;
    INTEGER(counts)[0] = 1;
    postold =  REAL(logmarg)[m] + log(REAL(priorprobs)[m]);
    memcpy(modelold, model, sizeof(int)*p);
    m = 0;
    int *varin= ivecalloc(p);
    int *varout= ivecalloc(p);
    double problocal = REAL(plocal)[0];
    


  Rprintf("Now Sample Models with MCMC  \n");  


   m = 0;

  while (nUnique < k && m < INTEGER(BURNIN_Iterations)[0]) {
    
    memcpy(model, modelold, sizeof(int)*p);
    pmodel =  n_sure;
    
    MH = GetNextModelCandidate(pmodel_old, n, n_sure, model, vars, problocal,
                               varin, varout, Rparents);
    
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
      prior_m = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
      if (prior_m == 0.0) {
        MH *= 0.0;
      }
      else {
        new_loc = nUnique;
        PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
        PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
        PROTECT(Rse_m = NEW_NUMERIC(pmodel));
        model_m = GetModel_m(Rmodel_m, model, p);
        
        R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY, pmodel, p, nobs, m, &mse_m,
                        &rank_m, pivot, tol);
        gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
        
        postnew = logmargy + log(prior_m);
        MH *= exp(postnew - postold);
      }}
    else {
      new_loc = branch->where;
      postnew =  REAL(logmarg)[new_loc] +
        log(REAL(priorprobs)[new_loc]);
      MH *=  exp(postnew - postold);
    }
    
    //    Rprintf("MH new %lf old %lf\n", postnew, postold);
    if (unif_rand() < MH) {
      if (newmodel == 1) {
        if ((m % thin) == 0 )  {
          
          new_loc = nUnique;
          insert_model_tree(tree, vars, n, model, nUnique);
          INTEGER(modeldim)[nUnique] = pmodel;
          INTEGER(rank)[nUnique] = rank_m;
          
          //record model data
          SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, nUnique);
          SetModel(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2,nUnique);
          
          ++nUnique;
         }
        else UNPROTECT(3);
       }
      
      old_loc = new_loc;
      postold = postnew;
      pmodel_old = pmodel;
      memcpy(modelold, model, sizeof(int)*p);
      
    } else  {
      if (newmodel == 1 && prior_m > 0) UNPROTECT(3);
    }
    
   if ( (m % thin) == 0) {
      
      INTEGER(counts)[old_loc] += 1;
      
      for (i = 0; i < n; i++) {
        // store in opposite order so nth variable is first
        real_model[n-1-i] = (double) modelold[vars[i].index];
        REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
      }
      nsamples++;
    }
    m++;
  }
  
  
  // Compute MCMC inclusion probabilities

  for (i = 0; i < n; i++) {
     REAL(MCMCprobs)[vars[i].index] /= (double) m;
  }
  Rprintf("\n Num Unique models %d  in %d MCMC iterations \n", nUnique, m);


// Compute marginal probabilities
  mcurrent = nUnique;
  compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
  compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);



//  Now sample W/O Replacement
 Rprintf("Sample w/out Replacement Now \n", nUnique);
 
 INTEGER(NumUnique)[0] = nUnique;

 m = nUnique; 
 
 if (nUnique < k) {
   update_probs(probs, vars, mcurrent, k, p);
   Rprintf("updating tree for SWOR\n");
   update_tree(modelspace, tree, modeldim, vars, k,p,n,mcurrent, modelwork);
   Rprintf("Done!\n");
   
  for (m = nUnique;  m < k && lessThanOne(pigamma[0]); m++) {
    INTEGER(modeldim)[m] = 0;
    for (i = n; i < p; i++)  {
      INTEGER(modeldim)[m]  +=  model[vars[i].index];
    }

    branch = tree;

    GetNextModel_swop(branch, vars, model, n, m, pigamma, problocal,
                      modeldim, bestmodel, Rparents);
    
    /* Now subtract off the visited probability mass. */
    branch=tree;
    Substract_visited_probability_mass(branch, vars, model, n, m, pigamma,eps);
    
    /* Now get model specific calculations */
    pmodel = INTEGER(modeldim)[m];
    PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
    memset(INTEGER(Rmodel_m), 0, pmodel * sizeof(int));
    PROTECT(Rcoef_m = NEW_NUMERIC(pmodel));
    PROTECT(Rse_m = NEW_NUMERIC(pmodel));
    model_m = GetModel_m(Rmodel_m, model, p);
    
    R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork, yty, SSY,
                    pmodel, p, nobs, m, &mse_m, &rank_m, pivot, tol);
    INTEGER(rank)[m] = rank_m;
    // initialize
    logmargy= 0.0;
    shrinkage_m = 1.0;
    gexpectations(p, rank_m, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
    //    Rprintf("rank %d dim %d\n", rank_m, pmodel);
    //		gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
    
    prior_m = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
    SetModel2(logmargy, shrinkage_m, prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
    SetModel_lm(Rcoef_m, Rse_m, Rmodel_m, mse_m, R2_m,	beta, se, modelspace, mse, R2,m);
    UNPROTECT(3);
    
    REAL(sampleprobs)[m] = pigamma[0];
    
/*
    if (m > 1) {
      rem = modf((double) m/(double) update, &mod);
      if (rem  == 0.0) {
      	mcurrent = m;
	      compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
	      compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);
	      if (update_probs(probs, vars, mcurrent, k, p) == 1) {
        	  Rprintf("Updating Model Tree %d \n", m);
	           update_tree(modelspace, tree, modeldim, vars, k, p, n ,mcurrent,
                         modelwork);
	           Rprintf("Done with update\n");
	        }
        }
      }
 */

  // end SWOR  
  }
 }
 
 /*. for when add heridity 
 if (m < k) {  // resize
   k = m;
   SETLENGTH(modelspace, m);
   SETLENGTH(logmarg, m);
   SETLENGTH(modelprobs, m);
   SETLENGTH(priorprobs, m);
   SETLENGTH(sampleprobs, m);
   SETLENGTH(beta, m);
   SETLENGTH(se, m);
   SETLENGTH(mse, m);
   SETLENGTH(shrinkage, m);
   SETLENGTH(modeldim, m);
   SETLENGTH(R2, m);
   SETLENGTH(rank, m);
   //  	Rprintf("m %d k %d", m, LENGTH(modelprobs));
 }
*/  
 

  Rprintf("Done with sampling - summaries m = %ld k = %ld \n", m, k);
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
  
  SET_VECTOR_ELT(ANS, 12, rank);
  SET_STRING_ELT(ANS_names, 12, mkChar("rank"));

  SET_VECTOR_ELT(ANS, 13, counts);
  SET_STRING_ELT(ANS_names, 13, mkChar("freq"));

  SET_VECTOR_ELT(ANS, 14, MCMCprobs);
  SET_STRING_ELT(ANS_names, 14, mkChar("probs.MCMC"));

  SET_VECTOR_ELT(ANS, 15, NumUnique);
  SET_STRING_ELT(ANS_names, 15, mkChar("n.Unique"));

  setAttrib(ANS, R_NamesSymbol, ANS_names);
 
  Rprintf("reset seed\n");
  PutRNGstate();
  
  Rprintf("free protected\n");
  UNPROTECT(nProtected);
  
  Rprintf("Return\n");
  return(ANS);
}




