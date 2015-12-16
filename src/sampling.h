#include <stdio.h>
#include <stdlib.h>
#include "mconf.h"
#include <math.h> 
#include <string.h>
#include <float.h>
#include <R.h> 
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


/* Defines. */
#define TRUE 1
#define FALSE 0
#define Bit unsigned char

/* Structs. */
struct Var {
  double prob;
  double logit;
  char flip;
  char leaveout;
  int index;
};

typedef struct Node *NODEPTR;

struct Node {
  double prob;
  int update;
  int counts_1;
  int counts_0;
  double logmarg;
  int where;
  NODEPTR zero;
  NODEPTR one;  
};



/* Subroutines. */
double hyp2f1(double a, double b, double c, double x);
void compute_margprobs(SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, double *margprobs, int k, int p);
void compute_margprobs_file(SEXP modeldim, SEXP Rmodelprobs, double *margprobs, int k, int p, FILE *file, int *model);
double beta_binomial(int modeldim, int p, double *hyper);
double Bernoulli(int *model, int p, double *hyper);
double compute_prior_probs(int *model, int modeldim, int p, SEXP modelprior);
void compute_margprobs_old(Bit **models, SEXP Rmodelprobs, double *margprobs, int k, int p);
void compute_modelprobs(SEXP modelprobs, SEXP logmarg, SEXP priorprobs,  int k);
void set_bits(char *bits, int subset, int *pattern, int *position, int n);
int compare(struct Var *i, struct Var *j);
	/* For sim. */
double *makeprob(double *prob, double *y, double c, double w, double
		 sigma2, int p, int inc_int);
void update_tree(SEXP modelspace, struct Node *tree, SEXP modeldim, struct Var *vars, int k, int p, int n, int kt, int *model);
void update_tree_file(struct Node *tree, SEXP modeldim, struct Var *vars, int k, int p, int n, int kt, FILE *file);
double random_switch(int *model, struct Var *vars, int n, int pmodel, int *varin, int *varout);
double random_walk(int *model, struct Var *vars, int n);
int update_probs(double *probs, struct Var *var, int m, int k, int p);
void update_MCMC_probs(double *probs, struct Var *vars, int n, int p);
void print_subset(int subset, int rank, Bit **models, Bit *model,
		  double *subsetsum, int *pattern, int *position, 
		  int n, struct Var *vars, int p);
int topk(Bit **models, double *prob, int k, struct Var *vars, int n, int p);
void do_insert(int child, double *subsetsum, int *queue);		
int get_next(double *subsetsum, int *queue, int *queuesize);
void insert_children(int subset, double *list, double *subsetsum, 
		     int *queue, int *queuesize, int *tablesize,
		     int *parent, int *pattern, int *position,
		     int *type, char *bits, int n);
void set_bits(char *bits, int subset, int *pattern, int *position, int n);
int sortvars(struct Var *vars, double *prob, int p);

void Lapack_chol2inv(double *cov, int p,double *covwork);
void F77_NAME(ch2inv)(double *cov, int *p, int *nr, double *covwork, int *info);
void F77_NAME(dpofa)(double *a, int *lda, int *n, int *info);
void F77_NAME(dposl)(double *a, int *lda, int *n, double  *b);
void F77_NAME(dpoco)(double *a, int *lda, int *n, double *rcond, double *z, int  *info);
void F77_NAME(dpodi)(double *a, int *lda, int *n, double *det,int  *job);
void F77_NAME(drsk)(char *UPLO, char *TRANS, int *N, int *K, double *ALPHA,
           double *A, int *LDA, double *BETA, double *C, int *LDC);

void  F77_NAME(dqrls)(double *x, int *n, int *p, double *y, int *inc, double *tol,
	    double *start, double *residuals, double *effects, int *rank, 
	    int *pivot, double *qraux, double *work);

//void F77_NAME(dgtsv)( int *n, int *nrhs, double *dl, double *d, double *du, double *b, double *ldb, int *info);

/*  Following defined in R_ext/BLAS.h 

double ddot_(int *n, double *x, int *incx, double *y, int *incy);
void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
void dsymv_(const char *uplo, int *n, double *alpha, double *A, int *lda, double *X, int *incx, double *beta, double *Y, int *incy);
void dsyrk_(const char *uplo, const char *trans,
		const int *n, const int *k,
		const double *alpha, const double *a, const int *lda,
		const double *beta, double *c, const int *ldc);
void dgemv_(const char *trans, int *m, int *n, double *alpha, double *A, 
            int *lda, double  *x, int *incx, double *beta, double *y,
            int *incy);

*/

void cholregold(double *residuals, double *X, double *XtX, double *coefficients,double *se, double *mse,  int p, int n);
void cholreg(double *XtY, double *XtX, double *coefficients,double *se, double *mse,  int p, int n);
double quadform (double *bwork, double *R,  int p);
double **matalloc(int das,int dbs);
int **imatalloc(int das,int dbs);
unsigned char **cmatalloc(int das,int dbs);
void  freechmat(unsigned char **mat, int  nr);
void  freemat(double **mat, int  nr);
double *vecalloc(int das);
int *ivecalloc(int das);
int withprob(double p);
void olsreg(double *Y, double *Xmat, double *coefficients, double *se, double *mse, int *p, int *n, int *pivot,double *qraux, double *work, double *residuals, double *effects, double *v, double *betaols);
void gexpectations(int p, int pmodel, int nobs, double R2, double alpha,
int method, double R2Full, double SSY, double *logmarg, double *shrinkage); 
double LogBF_Hg_null(double r2curr, int n, int d, double alpha, int gpower);
void posroot(double a, double b, double c, double *root, double *status);
double lik_null_HG(double g, double R2,int n,int k, double alpha, int gpower);
double info_null_HG(double g, double R2, int n, int k, double alpha);
double logBF_gprior(double Rsquare, int n,  int p, double g);
double logBF_hyperGprior(double R2, int n, int p, double alpha);
double logBF_hyperGprior_laplace(double R2,  int n, int p, double alpha);
double shrinkage_hyperg(double Rsquare, int n, int p, double alpha);
double shrinkage_laplace(double Rsquare, int n, int p, double alpha);
double shrinkage_EB_local(double Rsquare, int n, int p, double alpha);
double logBF_intrinsic_hyperGprior(double Rsquare, double RSquareFull, int n, int p, int pm, double alpha);
double logBF_EB(double Rsquare, int n, int p, double alpha);
double LogBF_ZS_full(double r2full, double r2curr, int n, int ptotal, int d);
double LogBF_ZS_null(double R2, int n, int d);
double E_ZS_approx_null(double R2, int n, int d);
double BIC(double Rsquare, int n,  int p, double SSY);
double AIC(double Rsquare, int n,  int p, double SSY);
NODEPTR make_node(double pr);
SEXP getListElement(SEXP list, char *str);
void PrecomputeData(double *Xwork, double *Ywork, double *wts, double **pXtXwork, double **pXtYwork, double **pXtX, double **pXtY, double *yty, double *SSY, int p, int nobs);

double CalculateRSquareFull(double *XtY, double *XtX, double *XtXwork, double *XtYwork, SEXP Rcoef_m, SEXP Rse_m, int p, int nobs, double yty, double SSY);
