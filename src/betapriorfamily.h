
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

extern double loghyperg1F1(double, double, double, int);
extern double shrinkage_chg(double a, double b, double Q, int laplace);
