#include "bas.h"

NODEPTR make_node(double pr) {
  NODEPTR newnode;
  newnode = (struct Node *) R_alloc(1, sizeof(struct Node));
  newnode->prob = pr;
  newnode->counts_1=0;
  newnode->counts_0=0;
  newnode->update = 0;
  newnode->logmarg = 0.0;
  newnode->where = -1;
  newnode->one = NULL;
  newnode->zero = NULL;
  return(newnode);
}

// not used now
// # nocov start
void deallocate_tree(struct Node *tree);

void deallocate_tree(struct Node *tree) {
  if (!tree) return;
  deallocate_tree(tree->one);
  deallocate_tree(tree->zero);
  Free(tree);
}
// # nocov end


void insert_model_tree(struct Node *tree, struct Var *vars,  int n, int *model, 
                       int num_models) {
  int i, bit;
  NODEPTR branch;

  branch = tree;

  // added bit to add counts //
  for (i = 0; i< n; i++) {
      bit =  model[vars[i].index];

      if (bit == 1) {
        branch->counts_1 += 1;
      	if (i < n-1 && branch->one == NULL)
	        branch->one = make_node(vars[i+1].prob);
      	if (i == n-1 && branch->one == NULL) {
	        branch->one = make_node(0.0);
      	}
      	branch = branch->one;
      }
      else {
        branch->counts_0 += 1;
	      if (i < n-1 && branch->zero == NULL)
	      branch->zero = make_node(vars[i+1].prob);
	      if (i == n-1 && branch->zero == NULL){
	      branch->zero = make_node(0.0);
	      }
	      branch = branch->zero;
      }
  }
  branch->where = num_models;
}

typedef int (*compfn)( const void* , const void*);

int sortvars(struct Var *vars, double *prob, int p)
{
	int i, n;

	/* Fill in variable information. */
	for (i = 0; i < p; i++) {
		vars[i].prob = prob[i];
		vars[i].index = i;
	}

	/* Make "list" from "probs".  Involves sorting and flipping and such. */
	n = 0;
	for (i = 0; i < p; i++) {
		if (vars[i].prob < 0.0) {
//			REprintf("Warning: Probability %d (%lf) less than zero, setting to zero.\n",
//				i, vars[i].prob);  
// # nocov start
			vars[i].leaveout = TRUE;
			vars[i].prob = 0.0;
// # nocov end			
		}
		else if (vars[i].prob == 0.0)
			vars[i].leaveout = TRUE;	/* Must be out. */
		else if (vars[i].prob < .5) {
			vars[i].leaveout = FALSE;
			vars[i].logit = log((1.0-vars[i].prob)/(vars[i].prob));
			vars[i].flip = TRUE;
			n++;
		}
		else if (vars[i].prob < 1.0) {
			vars[i].leaveout = FALSE;
			vars[i].logit = log((vars[i].prob)/(1.0-vars[i].prob));
			vars[i].flip = FALSE;
			n++;
		}
		else if (vars[i].prob == 1.0)
			vars[i].leaveout = TRUE;	/* Must be in. */
		else {
/*			REprintf("Warning: Probability %d (%lf) more than one, setting to one.\n",
				i, vars[i].prob); */
// # nocov start
			vars[i].leaveout = TRUE;
			vars[i].prob = 1.0;
// # nocov end			
		}
	}

/*	if (n == 0) {
		REprintf("Warning: prior inclusion probabilities are all 0 or 1 - Only 1 model!\n");
	}
 */
	/* Ok, vars is set up.  Need to sort to get "list". */
	qsort((char *) vars, p, sizeof(struct Var),(compfn) compare);


	return(n);
}


int compare(struct Var *i, struct Var *j)
{
  if (i->leaveout) return(1);
  if (j->leaveout) return(-1);
  if (i->logit > j->logit)
    return (-1);
  if (i->logit < j->logit)
    return (1);
  return (0);
}


int update_probs(double *probs, struct Var *vars, int m, int k, int p) {
	int i, update;
	double wt, newprob, diff;
	wt = (double) m / (double) k;
	update = 0;
	for (i=0, diff=0.0; i <p; i++){
		diff += (probs[vars[i].index] - vars[i].prob)*(probs[vars[i].index] - vars[i].prob);
	}
	if (sqrt(diff/ (double) p) >  .025) {
		update = 1;
		for (i = 0;   i < p; i++) {
			if (m < p) {
				newprob = probs[vars[i].index]* wt + vars[i].prob * (1 - wt);
			} else  newprob = probs[vars[i].index];
			if (newprob > .975) newprob = .975;
			if (newprob < .025) newprob = .025;
			vars[i].prob = newprob;
		}
	}
	return(update);
}


// used with glm_*.c
void SetModel1(SEXP Rfit, SEXP Rmodel_m,
	       SEXP beta, SEXP se, SEXP modelspace,
	       SEXP deviance, SEXP R2, SEXP Q, SEXP Rintercept, int m) {
	SET_ELEMENT(beta, m, getListElement(getListElement(Rfit, "fit"),"coefficients"));
	SET_ELEMENT(se, m, getListElement(getListElement(Rfit, "fit"),"se"));
	SET_ELEMENT(modelspace, m, Rmodel_m);

	REAL(R2)[m] = NA_REAL;
	REAL(deviance)[m] = REAL(getListElement(getListElement(Rfit, "fit"),"deviance"))[0];
	REAL(Q)[m] = REAL(getListElement(getListElement(Rfit, "lpy"),"Q"))[0];
	REAL(Rintercept)[m] = REAL(getListElement(getListElement(Rfit, "lpy"),"intercept"))[0];
}

void update_tree(SEXP modelspace, struct Node *tree, SEXP modeldim, 
                 struct Var *vars, int k, int p, int n, int kt, int *model)
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

void update_tree_AMC(SEXP modelspace, struct Node *tree, SEXP modeldim, 
                 struct Var *vars, int k, int p, int n, int kt, int *model, double *real_model, double *marg_probs, double *Cov, double delta)
{
  int i,m, bit;
  double prone, pigamma, przero;
  SEXP model_m;
  struct Node *branch;
  
  memset(model, 0, n*sizeof(int));
  memset(real_model, 0.0, n*sizeof(double));
  for (m = 0;  m <= kt; m++) {
    branch = tree;
    PROTECT(model_m = VECTOR_ELT(modelspace, m));
    
    for (i = 0; i < INTEGER(modeldim)[m]; i++)
      model[INTEGER(model_m)[i]] = 1;
    
    pigamma = 0.0;
    for (i = 0; i < n; i++) {
      real_model[i] = (double) model[vars[i].index];
      if (branch->update != kt) {
        branch->prob = cond_prob(real_model,i, n, marg_probs,Cov, delta);
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

void CreateTree_with_pigamma(NODEPTR branch, struct Var *vars,
                             int *bestmodel, int *model, int n,
                             int m, SEXP modeldim, double *pigamma,
                             SEXP Rparents) {
  double prob_parents;
  for (int i = 0; i< n; i++) {
    pigamma[i] = 1.0;
    int bit =  bestmodel[vars[i].index];
    model[vars[i].index] = bit;
    INTEGER(modeldim)[m]  += bit;
    if (bit == 1) {
      for (int j=0; j<=i; j++)  pigamma[j] *= branch->prob;
      if (i < n-1 && branch->one == NULL) {
        prob_parents = got_parents(bestmodel, Rparents, i+1, vars, n);
        branch->one = make_node(prob_parents);
      }
      if (i == n-1 && branch->one == NULL) {
        branch->one = make_node(0.0);
      }
      branch = branch->one;
    }
    else {
      for (int j=0; j<=i; j++)  pigamma[j] *= (1.0 - branch->prob);
      if (i < n-1 && branch->zero == NULL) {
        prob_parents = got_parents(bestmodel, Rparents, i+1, vars,n);
        branch->zero = make_node(prob_parents);
      }
      if (i == n-1 && branch->zero == NULL)
        branch->zero = make_node(0.0);
      branch = branch->zero;
    }
  }
}

void Substract_visited_probability_mass(NODEPTR branch, struct Var *vars, int *model, int n, int m,  double *pigamma, double eps) {
  for (int i = 0; i < n; i++) {
    int bit = model[vars[i].index];
    double prone = branch->prob;
    if (bit == 1) prone -= pigamma[i];
    double denom = 1.0 - pigamma[i];
    if (denom <= 0.0) {
/*      
      // # nocov start
      // should not be feasible
      if (denom < 0.0) {
        Rprintf("neg denominator %le %le %le !!!\n", pigamma, denom, prone);
        if (branch->prob < 0.0 && branch->prob < 1.0) {
          Rprintf("non extreme %le\n", branch->prob);
        }
        // # nocov end
      }
 */
      denom = 0.0;
    }
    else {
      if  (prone <= 0)  prone = 0.0;
      if  (prone > denom)  {
        if (prone <= eps) prone = 0.0;
        else prone = 1.0;
        /* Rprintf("prone > 1 %le %le %le %le !!!\n", pigamma, denom, prone, eps);*/
      }
      else prone = prone/denom;
    }
    // # nocov start
    // should not get here
    if (prone > 1.0 || prone < 0.0) {
      error("line 289: in tree-strutures.c sampling probability greater than 1\n");
    //  Rprintf("%d %d Probability > 1!!! %le %le  %le %le \n",
    //           m, i, prone, branch->prob, denom, pigamma);
    }
    // # nocov end
    branch->prob  = prone;
    if (bit == 1) branch = branch->one;
    else  branch = branch->zero;
  }
}

double GetNextModel_AMC(struct Var *vars,
                       int *model, int n, int m, SEXP modeldim,
                       SEXP Rparents, double *real_model, double*marg_probs, 
                       double *Cov, double delta) {
  double prob_parents = 1.0, pigamma = 1.0, prob = 0;
  int bit;
  
  for (int i = 0; i< n; i++) {
    
    prob = cond_prob(real_model,i, n, marg_probs,Cov, delta);
    bit = withprob(prob);
    model[vars[i].index] = bit;
    real_model[i] = (double) model[vars[i].index];
    
    if (bit == 1) {
        pigamma *= prob;  // calculate log probabilty of model
     }
    else {
      pigamma *= 1.0 - prob;
    }
   if (i < n-1) {
        //  check if parents
        prob_parents *=  got_parents(model, Rparents, i+1, vars, n);
      }
  }
  if (prob_parents <= 0) pigamma = 0;
  return(pigamma);
}

void GetNextModel_swop(NODEPTR branch, struct Var *vars,
                       int *model, int n, int m,  double *pigamma,
                       double problocal, SEXP modeldim, int *bestmodel,
                       SEXP Rparents) {
  double prob_parents = 1.0;
  for (int i = 0; i< n; i++) {
    pigamma[i] = 1.0;
    int bit =  withprob(branch->prob);
    
    model[vars[i].index] = bit;
    INTEGER(modeldim)[m]  += bit;
    
    if (bit == 1) {
      for (int j=0; j<=i; j++) {
        pigamma[j] *= branch->prob; } // calculate probabilty of model
      if (i < n-1 && branch->one == NULL) {
        //  add new branch
        //  check if parents
        prob_parents =  got_parents(model, Rparents, i+1, vars,n);
        branch->one = make_node(prob_parents);
      }
      if (i == n-1 && branch->one == NULL) branch->one = make_node(0.0);
      branch = branch->one;
    } else {
      for (int j=0; j<=i; j++)  pigamma[j] *= (1.0 - branch->prob);
      if (i < n-1 && branch->zero == NULL)
      {
        //  add new branch
        //  check if parents
        prob_parents =  got_parents(model, Rparents, i+1, vars, n);
        branch->zero = make_node(prob_parents);
      }
      if (i == n-1 && branch->zero == NULL) branch->zero = make_node(0.0);
      branch = branch->zero;
    }
  }
}

double got_parents(int *model, SEXP Rparents, int level, struct Var *var, int nsure)
{ double prob=1.0, *parents;
  int j=0, p, nsibs=0;
  
  int *dims = INTEGER(getAttrib(Rparents,R_DimSymbol));
  p = dims[0];
  if (p  == 1) {
    prob = var[level].prob;
  }
  else {
    parents = REAL(Rparents);
    // Rprintf("level %d\n", level);
    //
    // Check the variables that are always included first
    
    for (j = nsure, nsibs= 0, prob=1.0; j < p;  j++) {
      if ((parents[var[level].index + p*var[j].index]) == 1.0) {
        if (model[var[j].index] == 0) {
          // missing parent so probability of model is 0
          prob *= 0.0;  // # nocov
          }
        if (model[var[j].index] == 1) {
          // got parent so probability of variable is 1
          prob *= 1.0;
          nsibs += parents[var[j].index + p*var[level].index];
        }
      }
    }
    
    // now check the rest
    if (prob > 0.0) {
      //for (j=0, nsibs=0, prob=1.0; j < level; j++) {
      for (j=0; j < level; j++) {
        if ((parents[var[level].index + p*var[j].index]) == 1.0) {
          if (model[var[j].index] == 0) {
            // missing parent so probability of model is 0
            prob *= 0.0;}
          if (model[var[j].index] == 1) {
            // got parent so probability doesn't change
            prob *= 1.0;
            nsibs += parents[var[j].index + p*var[level].index];
          }
        }
        /*    Rprintf("%d pos %d, index %d, parents %lf, model %d nsibs %d, prob %lf\n",
         j, var[j].index,var[level].index,parents[var[level].index + p*var[j].index],
         model[var[j].index], nsibs, prob);
         */
      }
      if ((nsibs == 0) && (prob > 0.0))  prob = var[level].prob;
    }
  }
  
  
  //  Rprintf("Prob 1 = %lf\n", prob);
  return(prob);
}


int *GetModel_m(SEXP Rmodel_m, int *model, int p) {
  int *model_m = INTEGER(Rmodel_m);
  for (int j = 0, l=0; j < p; j++) {
    if (model[j] == 1) {
      model_m[l++] = j;
    }
  }
  return model_m;
}

void CreateTree(NODEPTR branch, struct Var *vars, int *bestmodel,
                int *model, int n, int m, SEXP modeldim, SEXP Rparents) {
  for (int i = 0; i< n; i++) {
    int bit =  bestmodel[vars[i].index];
    if (bit == 1) {
      if (i < n-1 && branch->one == NULL)
        branch->one = make_node(-1.0);
      if (i == n-1 && branch->one == NULL)
        branch->one = make_node(0.0);
      branch = branch->one;
    } else {
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
}
//caller is responsible to call UNPROTECT(2);
double FitModel(SEXP Rcoef_m, SEXP Rse_m, double *XtY, double *XtX, int *model_m,
                double *XtYwork, double *XtXwork, double yty, double SSY, int pmodel, int p,
                int nobs, int m, double *pmse_m, int *rank_m, int pivot, double tol) {
  
  double *coefficients = REAL(Rcoef_m);
  double *se_m = REAL(Rse_m);
  
  
  for (int j=0; j < pmodel; j++) { //subsetting matrix
    XtYwork[j] = XtY[model_m[j]];
    for  ( int i = 0; i < pmodel; i++) {
      XtXwork[j*pmodel + i] = XtX[model_m[j]*p + model_m[i]];
    }
  }
  *pmse_m = yty;
  memcpy(coefficients, XtYwork, sizeof(double)*pmodel);
  
  if (pivot == 1) {
    *rank_m = cholregpivot(XtYwork, XtXwork, coefficients, se_m, pmse_m, pmodel, nobs, pivot, tol);
  }
  else {
    cholreg(XtYwork, XtXwork, coefficients, se_m, pmse_m, pmodel, nobs);
    *rank_m = pmodel;
  }
  
  double R2_m = 1.0 - (*pmse_m * (double) ( nobs - *rank_m))/SSY;
  //	 Rprintf("R2 = %lf, mse = %lf, rank = %d\n", R2_m, *pmse_m, *rank_m);
  if (R2_m < 0.0 || *rank_m == 1) R2_m = 0.0;
  
  return R2_m;
}

void SetModel2(double logmargy, double shrinkage_m, double prior_m,
               SEXP sampleprobs, SEXP logmarg, SEXP shrinkage, SEXP priorprobs, int m) {
  REAL(sampleprobs)[m] = 1.0;
  REAL(logmarg)[m] = logmargy;
  REAL(shrinkage)[m] = shrinkage_m;
  REAL(priorprobs)[m] = prior_m;
}

void SetModel(SEXP Rcoef_m, SEXP Rse_m, SEXP Rmodel_m, double mse_m, double R2_m,
              SEXP beta, SEXP se, SEXP modelspace, SEXP mse, SEXP R2, int m) {
  
  SET_ELEMENT(beta, m, Rcoef_m);
  SET_ELEMENT(se, m, Rse_m);
  SET_ELEMENT(modelspace, m, Rmodel_m);
  
  REAL(R2)[m] = R2_m;
  REAL(mse)[m] = mse_m;
  
  UNPROTECT(3);
}

void SetModel_lm(SEXP Rcoef_m, SEXP Rse_m, SEXP Rmodel_m, double mse_m, double R2_m,
                 SEXP beta, SEXP se, SEXP modelspace, SEXP mse, SEXP R2, int m) {
  
  SET_ELEMENT(beta, m, Rcoef_m);
  SET_ELEMENT(se, m, Rse_m);
  SET_ELEMENT(modelspace, m, Rmodel_m);
  
  REAL(R2)[m] = R2_m;
  REAL(mse)[m] = mse_m;
  
}





