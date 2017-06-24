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

void deallocate_tree(struct Node *tree);

void deallocate_tree(struct Node *tree) {
  if (!tree) return;
  deallocate_tree(tree->one);
  deallocate_tree(tree->zero);
  Free(tree);
}



void insert_model_tree(struct Node *tree, struct Var *vars,  int n, int *model, int num_models) {
  int i, bit;
  NODEPTR branch;

  branch = tree;

  // added bit to add counts //
  for (i = 0; i< n; i++) {
      bit =  model[vars[i].index];

      if (bit == 1) {
        branch->counts_1 += 1;
	if (i < n-1 && branch->one == NULL)
	  branch->one = make_node(-1.0);
	if (i == n-1 && branch->one == NULL) {
	  branch->one = make_node(0.0);
	}
	branch = branch->one;
      }
      else {
        branch->counts_0 += 1;
	if (i < n-1 && branch->zero == NULL)
	  branch->zero = make_node(-1.0);
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
			REprintf("Warning: Probability %d (%lf) less than zero, setting to zero.\n",
				i, vars[i].prob);
			vars[i].leaveout = TRUE;
			vars[i].prob = 0.0;
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
			REprintf("Warning: Probability %d (%lf) more than one, setting to one.\n",
				i, vars[i].prob);
			vars[i].leaveout = TRUE;
			vars[i].prob = 1.0;
		}
	}

	if (n == 0) {
		error("Probabilities are all 0 or 1 - Quitting!\n");
	}
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


void update_MCMC_probs(double *probs, struct Var *vars, int n, int p) {
  int i;
  double  newprob;

    for (i = 0;   i < n; i++) {
      newprob = probs[vars[i].index];
      if (newprob > .975) newprob = .975;
      if (newprob < .025) newprob = .025;
      vars[i].prob = newprob;
    }
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
