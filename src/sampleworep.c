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

SEXP sampleworep(SEXP Y, SEXP X, SEXP Rweights, SEXP Rprobinit, SEXP Rmodeldim, SEXP incint, SEXP Ralpha,SEXP method, SEXP modelprior, SEXP Rupdate, SEXP Rbestmodel, SEXP Rbestmarg, SEXP plocal)
{
  SEXP   Rse_m, Rcoef_m, Rmodel_m; 

  SEXP   RXwork = PROTECT(duplicate(X)), RYwork = PROTECT(duplicate(Y));
  
  int  nProtected = 2;
  int  nModels=LENGTH(Rmodeldim);
  
  //  Rprintf("Allocating Space for %d Models\n", nModels) ;
  SEXP ANS = PROTECT(allocVector(VECSXP, 12)); ++nProtected;
  SEXP ANS_names = PROTECT(allocVector(STRSXP, 12)); ++nProtected;
  SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
  SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
  SEXP modeldim =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
  SEXP beta = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
  SEXP se = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
  SEXP mse = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP modelprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP priorprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP logmarg = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP sampleprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;

  double *Xwork, *Ywork, *wts, *coefficients,*probs, shrinkage_m, 
         SSY, yty, ybar, mse_m, *se_m, 
         R2_m, RSquareFull, alpha, prone, denom, logmargy;
  int nobs, p, k, i, j, m, n, l, pmodel, *xdims, *model_m, *bestmodel, local;
  int mcurrent,  update;
  double  mod, rem, problocal, *pigamma, eps, *hyper_parameters;
  double *XtX, *XtY, *XtXwork, *XtYwork;
  double one=1.0, zero=0.0;
 
  int inc=1, p2;
  int *model, bit, *modelwork;	
  char uplo[] = "U", trans[]="T";
  struct Var *vars;	/* Info about the model variables. */
  NODEPTR tree, branch;

  /* get dimsensions of all variables */


  nobs = LENGTH(Y);
  xdims = INTEGER(getAttrib(X,R_DimSymbol));
  p = xdims[1];
  k = LENGTH(modelprobs);
  update = INTEGER(Rupdate)[0];
  eps = DBL_EPSILON;
  problocal = REAL(plocal)[0];

  /* Extract prior on models  */
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));

  /*  Rprintf("n %d p %d \n", nobs, p);  */

  Ywork = REAL(RYwork);
  Xwork = REAL(RXwork);
  wts = REAL(Rweights);
 
 /* Allocate other variables.  */ 
  XtX  = (double *) R_alloc(p * p, sizeof(double));
  XtXwork  = (double *) R_alloc(p * p, sizeof(double));
  XtY = vecalloc(p);  
  XtYwork = vecalloc(p);

  /* initialize X matrix */
  for (j=0, l=0; j < p; j++) {
    for (i = 0; i < p; i++) {
      XtX[j*p + i] = 0.0;
    }
    for (i=0; i < nobs; i++) {
      //       Xmat[i][j] =  REAL(X)[l];
      //  Xwork[l] = Xmat[i][j];
       Xwork[l] = Xwork[l]*wts[i];
       l = l + 1; 
    } 
  }

  //  PROTECT(Rprobs = NEW_NUMERIC(p));
  //initprobs = REAL(Rprobinit);
 probs =  REAL(Rprobs);
  //memcpy(probs, initprobs,  p*sizeof(double));

 p2 = p*p;
 ybar = 0.0; SSY = 0.0; yty = 0.0; 

 
 F77_NAME(dsyrk)(uplo, trans, &p, &nobs, &one, &Xwork[0], &nobs, &zero, &XtX[0], &p); 


 for (i = 0; i< nobs; i++) {
   Ywork[i] = Ywork[i]*wts[i];
  }

 //  ybar = ybar/ (double) nobs;

  ybar = F77_NAME(ddot) (&nobs, &Ywork[0], &inc, &wts[0], &inc)/
         F77_NAME(ddot) (&nobs, &wts[0], &inc, &wts[0], &inc);
  yty = F77_NAME(ddot)(&nobs, &Ywork[0], &inc, &Ywork[0], &inc);
  SSY = yty - (double) nobs* ybar *ybar;

  
  F77_NAME(dgemv)(trans, &nobs, &p, &one, &Xwork[0], &nobs, &Ywork[0], &inc, &zero, &XtY[0],&inc);
  
  alpha = REAL(Ralpha)[0];

  vars = (struct Var *) R_alloc(p, sizeof(struct Var));
  n = sortvars(vars, probs, p); 

  /* Make space for the models and working variables. */ 

  /*  pivot = ivecalloc(p); 
  qraux = vecalloc(p);
  work =  vecalloc(2 * p);
  effects = vecalloc(nobs); 
  v =  vecalloc(p * p); 
  betaols = vecalloc(p);
  */

 
  pigamma = vecalloc(p);
  model = ivecalloc(p);
  modelwork= ivecalloc(p);


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
    RSquareFull = 1.0;
    if (nobs > p) {
    RSquareFull =  1.0 - (mse_m * (double) ( nobs - p))/SSY;
    }
    UNPROTECT(2);
  }

  /* fill in the sure things */
  for (i = n; i < p; i++)  {
      model[vars[i].index] = (int) vars[i].prob;
  }


  GetRNGstate();
  tree = make_node(vars[0].prob);


  //  Rprintf("For m=0, Initialize Tree with initial Model\n");  

  m = 0;
  bestmodel = INTEGER(Rbestmodel);
  for (i = n; i < p; i++)  {
      model[vars[i].index] = bestmodel[vars[i].index];
      INTEGER(modeldim)[m]  +=  bestmodel[vars[i].index];
   }
  /* Rprintf("Create Tree\n"); */
   branch = tree;

   for (i = 0; i< n; i++) {
      pigamma[i] = 1.0;
      bit =  bestmodel[vars[i].index];


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
   
    /* Now subtract off the visited probability mass. */
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
	  /* Rprintf("prone > 1 %le %le %le %le !!!\n", pigamma, denom, prone, eps);*/
	}
	else prone = prone/denom;
      }
      if (prone > 1.0 || prone < 0.0) 
	Rprintf("%d %d Probability > 1!!! %le %le  %le %le \n",
		m, i, prone, branch->prob, denom, pigamma);

      branch->prob  = prone;
      if (bit == 1) branch = branch->one;
      else  branch = branch->zero;

    }

    pmodel = INTEGER(modeldim)[m];
    PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
    model_m = INTEGER(Rmodel_m);


    /*    Rprintf("Now get model specific calculations \n"); */

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

    R2_m = 1.0;
    if (nobs > pmodel) {
      R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;
    }
    
    SET_ELEMENT(beta, m, Rcoef_m);
    SET_ELEMENT(se, m, Rse_m);

    REAL(R2)[m] = R2_m;
    REAL(mse)[m] = mse_m;

    gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull, SSY, &logmargy, &shrinkage_m);
    REAL(sampleprobs)[m] = 1.0;
    REAL(logmarg)[m] = logmargy;
    REAL(shrinkage)[m] = shrinkage_m;
    //    REAL(priorprobs)[m] = beta_binomial(pmodel,p,
    //    hyper_parameters);
    REAL(priorprobs)[m] = compute_prior_probs(model,pmodel,p, modelprior);
    REAL(Rbestmarg)[0] = REAL(logmarg)[m];
    UNPROTECT(3);

  /*   Rprintf("model %d max logmarg %lf\n", m, REAL(logmarg)[m]); */

    /*  Rprintf("Now Sample the Rest of the Models \n"); */

  for (m = 1;  m < k; m++) {
    for (i = n; i < p; i++)  {
      INTEGER(modeldim)[m]  +=  model[vars[i].index];
    }

    branch = tree;

    for (i = 0; i< n; i++) {
      pigamma[i] = 1.0;
      bit =  withprob(branch->prob);
      local = withprob(problocal);
      if ( local == 1 && (branch->prob < .999) && (branch->prob > 0.001))
	 bit = bestmodel[vars[i].index];
      else { 
	bit =  withprob(branch->prob); }

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
	  Rprintf("neg denominator %le %le %le !!!\n", pigamma, denom, prone);
	  if (branch->prob < 0.0 && branch->prob < 1.0)
	    Rprintf("non extreme %le\n", branch->prob);}
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

      /*      Rprintf("%d %d \n",  branch->done, n - i); */
      /*      if (log((double) branch->done) < (n - i)*log(2.0)) {
	if (bit == 1) branch = branch->one;
	else  branch = branch->zero;
      }
      else {
	    branch->one = NULL;
	    branch->zero = NULL; 
	    break; } */
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

  
/*    olsreg(Ywork, Xwork, coefficients, se_m, &mse_m, &pmodel, &nobs, pivot,qraux,work,residuals,effects,v,betaols);   */
    R2_m = 1.0;
    if (nobs > pmodel) {
    R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;
    }
    
    SET_ELEMENT(beta, m, Rcoef_m);
    SET_ELEMENT(se, m, Rse_m);

    REAL(R2)[m] = R2_m;
    REAL(mse)[m] = mse_m;

   gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0],  RSquareFull, SSY, &logmargy, &shrinkage_m);
   REAL(logmarg)[m] = logmargy;
   REAL(shrinkage)[m] = shrinkage_m;
   REAL(priorprobs)[m] = compute_prior_probs(model,pmodel,p, modelprior);
   //   REAL(priorprobs)[m] = beta_binomial(pmodel,p, hyper_parameters);
   if (REAL(logmarg)[m] > REAL(Rbestmarg)[0]) {
      for (i=0; i < p; i++) {
	bestmodel[i] = model[i];}
      REAL(Rbestmarg)[0] = REAL(logmarg)[m];
    }
    
       

    if (m > 1) {
      rem = modf((double) m/(double) update, &mod);
      if (rem  == 0.0) {
	mcurrent = m;
	compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);        
	if (update_probs(probs, vars, mcurrent, k, p) == 1) {
	  //	  Rprintf("Updating Model Tree %d \n", m);
	  update_tree(modelspace, tree, modeldim, vars, k,p,n,mcurrent, modelwork);     
	}

      }}
    UNPROTECT(3);  
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

  setAttrib(ANS, R_NamesSymbol, ANS_names);
  UNPROTECT(nProtected);

  PutRNGstate();

  return(ANS);  
}


void update_tree(SEXP modelspace, struct Node *tree, SEXP modeldim, struct Var *vars, int k, int p, int n, int kt, int *model)
{
	int i,m, bit;
	double prone, pigamma, przero;
	SEXP model_m;
	struct Node *branch;

	for (m = 0;  m <= kt; m++) {
		branch = tree;
		PROTECT(model_m = VECTOR_ELT(modelspace, m));
		for (i = 0; i < p; i++)  model[i] = 0;
		for (i = 0; i < INTEGER(modeldim)[m]; i++) 
			model[INTEGER(model_m)[i]] = 1;

		pigamma = 0.0;
		for (i = 0; i < n; i++) {
  			if (branch->update != kt) {
				branch->prob = vars[i].prob;
				branch->update = kt;
			}
  			bit = model[vars[i].index];
			if (bit ==  1)  {
				pigamma += log(branch->prob);
				branch = branch->one;
			} else {
				pigamma += log(1.0 - branch->prob);
				branch = branch->zero;
			}
		}

		branch = tree;
		for (i = 0; i < n; i++) {
			bit = model[vars[i].index];
			if (bit == 1) {
				prone = (branch->prob - exp(pigamma));
				przero = 1.0 - branch->prob;
				pigamma -= log(branch->prob);
			} else {
				prone = branch->prob;
				przero = 1.0 - branch->prob  - exp(pigamma);
				pigamma -= log(1.0 - branch->prob);
			}
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





