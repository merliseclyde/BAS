#include "bas.h"

// [[register]]
SEXP deterministic(SEXP Y, SEXP X, SEXP Rweights, SEXP Rprobinit,
                   SEXP Rmodeldim, SEXP incint, SEXP Ralpha,
                   SEXP method, SEXP modelprior, SEXP Rpivot)
{
  SEXP   RXwork = PROTECT(duplicate(X)), RYwork = PROTECT(duplicate(Y));
  int nProtected = 2;

  int  nModels=LENGTH(Rmodeldim);
  int  pivot = LOGICAL(Rpivot)[0];

  SEXP ANS = PROTECT(allocVector(VECSXP, 12)); ++nProtected;
  SEXP ANS_names = PROTECT(allocVector(STRSXP, 12)); ++nProtected;
  SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
  SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
  SEXP modeldim =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
  SEXP rank =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
  SEXP beta = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
  SEXP se = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
  SEXP mse = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP modelprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP priorprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP logmarg = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
  SEXP sampleprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;

  SEXP Rse_m, Rcoef_m, Rmodel_m;
  double *Xwork, *Ywork, *wts, *coefficients,*probs,
    SSY, yty, mse_m, *se_m, pigamma,
    R2_m, RSquareFull, alpha, logmarg_m, shrinkage_m;
  double *XtX, *XtY, *XtXwork, *XtYwork;
  int nobs, p, k, i, j, m, n, l, pmodel, *xdims, *model_m, *model;
  Bit **models;
  int rank_m;
  struct Var *vars;	/* Info about the model variables. */


  /* get dimsensions of all variables */

  nobs = LENGTH(Y);
  xdims = INTEGER(getAttrib(X,R_DimSymbol));
  p = xdims[1];
  k = LENGTH(modelprobs);

  Ywork = REAL(RYwork);
  Xwork = REAL(RXwork);
  wts = REAL(Rweights);


	PrecomputeData(Xwork, Ywork, wts, &XtXwork, &XtYwork, &XtX, &XtY, &yty, &SSY, p, nobs);


	alpha = REAL(Ralpha)[0];

  vars = (struct Var *) R_alloc(p, sizeof(struct Var));
  probs =  REAL(Rprobs);
  n = sortvars(vars, probs, p);


  /* Make space for the models and working variables. */


  models = cmatalloc(k,p);
  model = (int *) R_alloc(p, sizeof(int));

  k = topk(models, probs, k, vars, n, p);



  /* Fit Full model */
  if (nobs <= p) {RSquareFull = 1.0;}
  else {


    Rcoef_m = NEW_NUMERIC(p); PROTECT(Rcoef_m);
    Rse_m = NEW_NUMERIC(p); PROTECT(Rse_m);
    coefficients = REAL(Rcoef_m);  se_m = REAL(Rse_m);

    memcpy(coefficients, XtY,  p*sizeof(double));
    memcpy(XtXwork, XtX, p*p*sizeof(double));
    memcpy(XtYwork, XtY,  p*sizeof(double));

    mse_m = yty;

    rank_m =  cholregpivot(XtYwork, XtXwork, coefficients, se_m, &mse_m, p, nobs);

    RSquareFull =  1.0 - (mse_m * (double) ( nobs - rank_m))/SSY;
    UNPROTECT(2);
  }

  /* now fit all top k models */

  for (m=0; m < k; m++) {
      pmodel = 0;
      pigamma = 1.0;
      for (j = 0; j < p; j++) {
          model[j] = (int) models[m][j];
          pmodel += (int) models[m][j];
          pigamma *= (double)((int) models[m][j])*probs[j] +
	  (1.0 - (double)((int) models[m][j]))*(1.0 -  probs[j]);
      }

      REAL(sampleprobs)[m] = pigamma;
      INTEGER(modeldim)[m] = pmodel;
      Rmodel_m = NEW_INTEGER(pmodel); PROTECT(Rmodel_m);
      model_m = INTEGER(Rmodel_m);


      for (j = 0, l=0; j < p; j++) {
	if (models[m][j]) {
           model_m[l] = j;
           l +=1;  }
      }

      INTEGER(modeldim)[m] = pmodel;
      SET_ELEMENT(modelspace, m, Rmodel_m);
      UNPROTECT(1);

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
//      cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, pmodel, nobs);

//      R2_m = 1.0 - (mse_m * (double) ( nobs - pmodel))/SSY;

      R2_m = FitModel(Rcoef_m, Rse_m, XtY, XtX, model_m, XtYwork, XtXwork,
                      yty, SSY, pmodel, p, nobs, m, &mse_m, &rank_m,
                      pivot);
      INTEGER(rank)[m] = rank_m;

      SET_ELEMENT(beta, m, Rcoef_m);
      SET_ELEMENT(se, m, Rse_m);

      REAL(R2)[m] = R2_m;
      REAL(mse)[m] = mse_m;

      gexpectations(p, pmodel, nobs, R2_m, alpha, INTEGER(method)[0], RSquareFull,
                    SSY, &logmarg_m, &shrinkage_m);
      REAL(logmarg)[m] = logmarg_m;
      REAL(priorprobs)[m] = compute_prior_probs( model, pmodel,p, modelprior);
      REAL(shrinkage)[m] = shrinkage_m;
      UNPROTECT(2);
  }

  compute_modelprobs(modelprobs, logmarg, priorprobs, k);
  compute_margprobs_old(models, modelprobs, probs, k, p);

    /*    freechmat(models,k); */
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

  return(ANS);

}

/* Find top k model function for deterministic search. */

int topk(Bit **models, double *prob, int k, struct Var *vars, int n, int p)
{
  int rank; 	/* index */

  int tablesize;	/* Number of entries in the subset table. */
  int i, current, qsize, queuesize;
  Bit *model;		/* output for current model  */
  double *list;		/* The list of numbers (logits). */
  double *subsetsum;	/* Sum of some subset. */
  int *type;	/* How does subset differ from parent? */
  int *queue;	/* Subset number of item in priority queue. */
  int *parent;	/* Parent of given subset. */
  int *position;	/* Bit position of rightmost 1 bit. */
  int *pattern;	/* Subset that differs only in bit position. */
  char *bits;


  /* Ask for too many models? */
  if (log((double) k)/log((double) 2.0) > (double) n) {
    REprintf("Warning: Asked for %d models, but there are only 2^%d with non-zero prob.\n",
	k, n);
    k = 1<<n;
  }

  /*  Rprintf("Deterministic sample of (approx) top %d models\n", k); */

  list =  (double *) R_alloc(n, sizeof(double));
  for (i = 0; i < n; i++) list[i] = vars[i].logit;

  model = (unsigned char *) R_alloc(n, sizeof(Bit));

  qsize = 2*k;  /* Largest number of items in priority queue. */
  subsetsum = (double *) R_alloc(qsize, sizeof(double));
  parent = (int *) R_alloc(qsize, sizeof(int));
  type =(int *) R_alloc(qsize, sizeof(int));
  position = (int *) R_alloc(qsize, sizeof(int));
  pattern = (int *) R_alloc(qsize, sizeof(int));
  queue = (int *) R_alloc(qsize, sizeof(int));
  bits = (char *) R_alloc(n, sizeof(char));

  queuesize = 0;
  /* Store relevant information for root node. */
  for(i = 0; i < n; subsetsum[0] += list[i], i++)
    ;
   parent[0] = -1;
   position[0] = -1;
   type[0] = 0;
   tablesize = 1;
   pattern[0] = -1;

  print_subset(0, 0, models, model, subsetsum, pattern, position, n, vars, p);
  insert_children(0, list, subsetsum, queue, &queuesize, &tablesize,
		  parent,pattern, position,type, bits, n);

  for (rank = 1; rank < k-1; rank++) {
    current = get_next(subsetsum, queue, &queuesize);
    print_subset(current, rank, models, model, subsetsum, pattern,
		 position, n, vars, p);
    insert_children(current, list, subsetsum, queue, &queuesize, &tablesize,
		    parent, pattern, position, type, bits,n);
  }
  if (k > 1) print_subset(queue[0], k-1, models, model, subsetsum,
			  pattern, position, n, vars, p);	/* Last one. */

  return(k);
}

/* Put subset's appropriate children into priority queue.  The current
   version inserts no more than 2 entries (although it takes O(n)) to
   figure out which two). */

void insert_children(int subset, double *list, double *subsetsum,
		     int *queue, int *queuesize, int *tablesize,
		     int *parent, int *pattern, int *position,
		     int *type, char *bits, int  n)
{
  int current;

  set_bits(bits, subset, pattern, position, n);

  /* Rightmost bit is off. */
  if (bits[n-1] == 0) {
    /* Put the kid on the table.  (I mean "in".) */
    *tablesize = 1 + *tablesize;
    current = *tablesize;
    subsetsum[current] = subsetsum[subset]-list[n-1];
    type[current] = 1;
    position[current] = n-1;
    parent[current] = subset;
    pattern[current] = subset;

    /* Put the subset into the priority queue. */
    queue[*queuesize] = current;
    do_insert(*queuesize, subsetsum, queue);	/* Fix up heap. */
    *queuesize = 1+ *queuesize;
  }

  /* Where's the rightmost that is *on*? */
  if (position[subset] <= 0) return;	/* Can't shift it---too far over. */
  if (bits[position[subset]-1] == 1) return; /* No place to shift to. */

  /* Put the kid on the table.  (I mean "in".) */
  *tablesize = 1 + *tablesize;
  current = *tablesize;
  subsetsum[current] =
    subsetsum[subset]+list[position[subset]]-list[position[subset]-1];

  type[current] = 2;
  position[current] = position[subset]-1;
  parent[current] = subset;
  pattern[current] = pattern[subset];

  /* Put the subset into the priority queue. */
  queue[*queuesize] = current;
  do_insert(*queuesize, subsetsum, queue);	/* Fix up heap. */
  *queuesize = 1+ *queuesize;

  return;

}

/* Fix up heap after insert. */
void do_insert(int child, double *subsetsum, int *queue)
{
  int swap, parent;

  /* Done if we get to the top. */
  while (child != 0) {
    parent = (child-1)/2;
    if (subsetsum[queue[child]] < subsetsum[queue[parent]]) return;

    /* Otherwise, things are out of order and need to be swapped. */
    swap = queue[child];
    queue[child] = queue[parent];
    queue[parent] = swap;

    child = parent;
  }
}

/* Get next subset from queue.  Returns top of queue.  Deletes
   top element from queue. */
int get_next(double *subsetsum, int *queue, int *queuesize)
{
  int parent, big_kid, kid1, kid2, swap, current;

  current = queue[0];
  queue[0] = queue[*queuesize-1];
  queue[*queuesize-1] = 0;	/* Just for neatness. */
  *queuesize = *queuesize -1;

  for (parent = 0; parent < *queuesize;) {
    kid1 = parent*2+1;
    kid2 = kid1+1;

    if (kid1 >= *queuesize) return(current);	/* At a leaf. */

    if (kid2 >= *queuesize) { 	/* Only one child. */
      if (subsetsum[queue[parent]] < subsetsum[queue[kid1]]) big_kid = kid1;
      else return(current);		/* Parent is best. */
    } else {			/* Two kids. */
      if ((subsetsum[queue[parent]] > subsetsum[queue[kid1]])&&
	  (subsetsum[queue[parent]] > subsetsum[queue[kid2]]))
	return(current);			/* Parent is bigger. */
      else {
	if (subsetsum[queue[kid1]] <= subsetsum[queue[kid2]]) big_kid = kid2;
	else big_kid = kid1;
      }
    }
    /* Now big_kid is the one to swap with. */
    swap = queue[parent];
    queue[parent] = queue[big_kid];
    queue[big_kid] = swap;

    parent = big_kid;
  }
  return(current);
}

/* Given an array of bits, and a subset, set the bits as approriate. */
void set_bits(char *bits, int subset, int *pattern, int *position, int n)
{int i;

  /* We know the bit pattern for ROOT. */
  for (i = 0; i < n; i++) bits[i] = 0;

  for (;subset != 0; subset = pattern[subset])
    bits[position[subset]] = 1;
}

/* Print current subset. */
void print_subset(int subset, int rank, Bit **models, Bit *model,
		  double *subsetsum, int *pattern, int *position,
		  int n, struct Var *vars, int p)
{
  int i;

  for (i = 0; i < n;  model[i] = 1, i++);
  for (i = subset; i!= 0; i = pattern[i])
    model[position[i]] = 0;

  for (i = 0; i < p; i++) {
    if (vars[i].leaveout)
      models[rank][vars[i].index] = (int) vars[i].prob;
    else if (vars[i].flip)
      models[rank][vars[i].index] = 1-model[i];
    else
      models[rank][vars[i].index] = model[i];
  }
}


int withprob(double p)
{
  if (unif_rand()<=p) return 1; else return 0;
}
