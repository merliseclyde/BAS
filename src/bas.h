// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
// before any R headers, or define in PKG_CPPFLAGS
#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mconf.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Constants.h>
#include <R_ext/Applic.h>
#include <R_ext/Rdynload.h>

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

double CalculateRSquareFull(double *XtY, double *XtX, double *XtXwork, double *XtYwork,
                            SEXP Rcoef_m, SEXP Rse_m, int p, int nobs, double yty, double SSY);
int *GetModel_m(SEXP Rmodel_m, int *model, int p);
void SetModel2(double logmargy, double shrinkage_m, double prior_m,
               SEXP sampleprobs, SEXP logmarg, SEXP shrinkage, SEXP priorprobs, int m);
void SetModel(SEXP Rcoef_m, SEXP Rse_m, SEXP Rmodel_m, double mse_m, double R2_m,
              SEXP beta, SEXP se, SEXP modelspace, SEXP mse, SEXP R2, int m);
void SetModel_lm(SEXP Rcoef_m, SEXP Rse_m, SEXP Rmodel_m, double mse_m, double R2_m,
                 SEXP beta, SEXP se, SEXP modelspace, SEXP mse, SEXP R2, int m);

double hyp2f1(double a, double b, double c, double x);
void compute_margprobs(SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, double *margprobs, int k, int p);
void compute_margprobs_file(SEXP modeldim, SEXP Rmodelprobs, double *margprobs, int k, int p, FILE *file, int *model);
double beta_binomial(int modeldim, int p, double *hyper);
double trunc_beta_binomial(int modeldim, int p, double *hyper);
double trunc_poisson(int modeldim, int p, double *hyper);
double trunc_power_prior(int modeldim, int p, double *hyper);
double Bernoulli(int *model, int p, double *hyper);
int no_prior_inclusion_is_1(int p, double *probs);
double compute_prior_probs(int *model, int modeldim, int p, SEXP modelprior, int noInclusionIs1);
void compute_margprobs_old(Bit **models, SEXP Rmodelprobs, double *margprobs, int k, int p);
void compute_modelprobs(SEXP modelprobs, SEXP logmarg, SEXP priorprobs,  int k);
void compute_modelprobs_HT(SEXP Rmodelprobs,  SEXP Rlogmarg, SEXP Rpriorprobs, SEXP Rsampleprobs, 
                           int k, int MC);
void set_bits(char *bits, int subset, int *pattern, int *position, int n);
int compare(struct Var *i, struct Var *j);
	/* For sim. */
double *makeprob(double *prob, double *y, double c, double w, double
		 sigma2, int p, int inc_int);
void update_tree(SEXP modelspace, struct Node *tree, SEXP modeldim, 
                 struct Var *vars, int k, int p, int n, int kt, int *model);
void update_tree_AMC(SEXP modelspace, struct Node *tree, SEXP modeldim, 
                 struct Var *vars, int k, int p, int n, int kt, int *model, 
                 double *real_model, double *marg_probs, double *Cov, double delta);
  
void update_tree_file(struct Node *tree, SEXP modeldim, struct Var *vars, int k, int p, int n, int kt, FILE *file);
double random_switch_heredity(int *model, struct Var *vars, int n, int pmodel, int *varin, int *varout, SEXP Rparents);
double random_walk_heredity(int *model, struct Var *vars, int n, SEXP Rparents);
double random_switch(int *model, struct Var *vars, int n, int pmodel, int *varin, int *varout);
double random_walk(int *model, struct Var *vars, int n);
int update_probs(double *probs, struct Var *var, int m, int k, int p);
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

int cholregpivot(double *XtY, double *XtX, double *coefficients,double *se, double *mse,  int p, int n, int pivot, double tol);
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
double ZS_logmarg(double R2, int n, int d, double rscale);
double ZS_shrinkage(double R2, int n, int d, double rscale);
void ZS_density(double *x, int n, SEXP Rtheta);
void ZS_density_shrinkage(double *x, int n, SEXP Rtheta);
double find_mode_g_JZS(double R2, int n, int d, double *root, double *status);
double E_ZS_approx_null(double R2, int n, int d);
void phi1(double *a, double *b, double *c, double *x, double *y, int *div, double *scale, double*phi, 
             int *npara);
void phi1_density(double *u, int n, SEXP Rtheta);
double tcch_int(double a, double b, double r, double s, double v,  double k);
void tcch_density(double *u, int n, SEXP Rtheta);
void tcch(double *a, double *b, double *r, double *s, double *v, double *theta, double *tcch, int *npara);
double phi1_int(double a, double b, double c, double x, double y, int div, double scale);
double BIC(double Rsquare, int n,  int p, double SSY);
double AIC(double Rsquare, int n,  int p, double SSY);
NODEPTR make_node(double pr);
SEXP getListElement(SEXP list, char *str);
void PrecomputeData(double *Xwork, double *Ywork, double *wts, double **pXtXwork, double **pXtYwork, double **pXtX, double **pXtY, double *yty, double *SSY, int p, int nobs);

double CalculateRSquareFull(double *XtY, double *XtX, double *XtXwork, double *XtYwork, SEXP Rcoef_m, SEXP Rse_m, int p, int nobs, double yty, double SSY);


//struct glmsfamily * make_glm_family(SEXP family);

typedef struct betapriorfamilystruc {
  const char *priorfamily;
  const char *samplingmodel;
  const char *priorclass;
  SEXP hyperparams;
  double (*logmarglik_fun)(SEXP hyperparams, int pmodel, double W, double loglik_mle, double logdet_Iintercept, int Laplace);
  double (*shrinkage_fun)(SEXP hyperparams, int pmodel, double W, int Laplace);
} betapriorptr;

betapriorptr * make_betaprior_structure(SEXP betaprior, SEXP glmfamily);

// CCH family
double CCH_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglike_mle, double logdet_Iintercept, int Laplace);
double CCH_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace);
// Robust

double robust_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglik_mle, double logdet_Iintercept, int Laplace );
double robust_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace );

double tCCH_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglik_mle, double logdet_Iintercept, int Laplace );
double tCCH_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace );
double HyperTwo(double a, double b, double c, double x, double y);

double intrinsic_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglik_mle, double logdet_Iintercept, int Laplace );
double intrinsic_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace );
double Jeffreys_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglike_mle, double logdet_Iintercept, int Laplace);
double TG_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglik_mle, double logdet_Iintercept, int Laplace );
double TG_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace );

double EB_local_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglik_mle, double logdet_Iintercept, int Laplace );
double EB_local_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace );

double betaprime_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglik_mle, double logdet_Iintercept, int Laplace );
double betaprime_glm_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace );

// BIC, AIC,...
double IC_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglike_mle, double logdet_Iintercept, int Laplace);
double IC_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace);
// g-prior
double g_prior_glm_logmarg(SEXP hyperparams, int pmodel, double W, double loglike_mle, double logdet_Iintercept, int Laplace);
double g_prior_shrinkage(SEXP hyperparams, int pmodel, double W, int Laplace);
// test-based BF
double testBF_prior_glm_logmarg(SEXP hyperparams, int pmodel, double W,
                                double loglik_mle, double logdet_Iintercept,
                                int Laplace );

extern double loghyperg1F1(double, double, double, int);
extern double shrinkage_chg(double a, double b, double Q, int laplace);

#ifndef R_STATS_FAMILY_H
#define R_STATS_FAMILY_H

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif

//struct glmsfamily * make_glm_family(SEXP family);

typedef struct glmfamilystruc {
  const char *family;
  const char *link;
  void (*mu_eta)(double *eta, double *mu, int n);
  void (*linkfun)(double *mu, double *eta, int n);
  void (*variance)(double * mu, double *var, int n);
  void (*dev_resids)(double *y, double *mu, double *weights, double *resids, int n);
  void (*linkinv)(double *eta, double *mu, int n);
  void (*initialize)(double *Y, double *mu, double *weights, int n);
  double (*dispersion)(double *resid,  double *weights, int n, int rank);
  void (*info_matrix) (double *y, double *mu, double *weights, double *var, int n);
  double (*loglik) (double *y, double *mu, double *weights, double devb, int n);
} glmstptr;

glmstptr * make_glmfamily_structure(SEXP family);

// Binomial family
double binomial_loglik(double *Y, double *mu, double *wts, double devb, int n);
void logit_link(double *mu, double *eta, int n);
void logit_variance(double *mu, double *var, int n);
void logit_linkinv(double *eta, double *mu, int n);
void logit_mu_eta(double *eta, double *mu, int n);
void logit_info(double *y, double *mu, double *weights, double *var, int n);
void binomial_dev_resids(double *y, double *mu, double *wt, double *res, int n);
double binomial_dispersion(double *resid,  double *weights, int n, int rank);
void binomial_initialize(double *Y, double *mu,  double *weights, int n);
// Poisson
double poisson_loglik(double *Y, double *mu, double *wts, double devb, int n);
double poisson_dispersion(double *resid,  double *weights, int n, int rank) ;
void log_link(double *mu, double *eta, int n);
void poisson_variance(double *mu, double *var, int n);
void log_linkinv(double *eta, double *mu, int n);
void log_mu_eta(double *eta, double *mu, int n);
void poisson_log_info(double *y, double *mu, double *weights, double *var, int n);
void poisson_dev_resids(double *y, double *mu, double *wt, double *res, int n);
void poisson_initialize(double *Y, double *mu,  double *weights, int n);
// Gamma
double gamma_loglik(double *Y, double *mu, double *wts, double devb, int n);
void gamma_variance(double *mu, double *var, int n);
void gamma_dev_resids(double *y, double *mu, double *wt, double *res, int n);
void gamma_initialize(double *Y, double *mu,  double *weights, int n);


// Normal
double Gaussian_dispersion(double *resid,  double *weights, int n, int rank) ;



// General Functions
double deviance(double *res, int n);
void chol2se(double *qr, double *se, double *cov, double *covwork, int p, int n);
void QR2cov(double *qr, double *cov, double *covwork, int p, int n);
#endif

void  update_Cov(double *Cov, double *priorCov, double *SSgam, double *marg_probs, double lambda, int n, int m, int print);
double cond_prob(double *model, int j, int n, double *mean, double *beta_matrix , double delta) ;

void insert_model_tree(struct Node *tree, struct Var *vars,  int n, int *model, int num_models);

int *GetModel_m(SEXP Rmodel_m, int *model, int p);

void CreateTree(NODEPTR branch, struct Var *vars, int *bestmodel,
                int *model, int n, int m, SEXP modeldim, SEXP Rparents);

void CreateTree_with_pigamma(NODEPTR branch, struct Var *vars, int *bestmodel, int *model, int n, int m,
                             SEXP modeldim, double *pigamma, SEXP Rparents);

double GetNextModelCandidate(int pmodel_old, int n, int n_sure, int *model, struct Var *vars, double problocal,
                             int *varin, int *varout, SEXP Rparents);
double got_parents(int *model, SEXP Rparents, int level, struct Var *var, int nsure);

void GetNextModel_swop(NODEPTR branch, struct Var *vars, int *model, int n, int m,  double *pigamma,
                       double problocal, SEXP modeldim,int *bestmodel,
                       SEXP Rparents);
double GetNextModel_AMC(struct Var *vars,
                      int *model, int n, int m, SEXP modeldim,
                      SEXP Rparents, double *real_model, double*marg_probs, 
                      double *Cov, double delta);
void Substract_visited_probability_mass(NODEPTR branch, struct Var *vars, int *model, int n, int m, double *pigamma, double eps);

void SetModel1(SEXP Rfit, SEXP Rmodel_m,
               SEXP beta, SEXP se, SEXP modelspace, SEXP deviance, SEXP R2, SEXP Q, SEXP Rintercept, int m);

void SetModel2(double logmargy, double shrinkage_m, double prior_m,
               SEXP sampleprobs, SEXP logmarg, SEXP shrinkage, SEXP priorprobs, int m);
double FitModel(SEXP Rcoef_m, SEXP Rse_m, double *XtY, double *XtX, int *model_m,
                double *XtYwork, double *XtXwork, double yty, double SSY, int pmodel, int p,
                int nobs, int m, double *pmse_m, int *rank_m, int pivot, double tol);
SEXP glm_FitModel(SEXP RX, SEXP RY, SEXP Rmodel_m, //input data
                  SEXP Roffset, SEXP Rweights,
                  glmstptr * glmfamily, SEXP Rcontrol,
                  SEXP Rlaplace, betapriorptr * betapriorfamily);

SEXP glm_bas(SEXP RX, SEXP RY, glmstptr * family, SEXP Roffset, SEXP Rweights, SEXP Rcontrol);

SEXP gglm_lpy(SEXP RX, SEXP RY,SEXP Rcoef, SEXP Rmu, SEXP Rdeviance, SEXP Rweights, glmstptr * glmfamily, betapriorptr * betapriorfamily, SEXP Rlaplace);


// issue 38
static inline int lessThanOne(double a)
{
  // DBL_EPSILON is too restrictive. This might need further tweaking
  double LOCAL_DBL_EPSILON = 1E-10;
  return (1.0 - a) >= (LOCAL_DBL_EPSILON);
}

