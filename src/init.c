#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>



/* .C calls */
extern void gexpectations_vect(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hypergeometric1F1(void *, void *, void *, void *, void *, void *);
extern void hypergeometric2F1(void *, void *, void *, void *, void *);
extern void logHyperGauss2F1(void *, void *, void *, void *, void *);
extern void phi1(void *, void *, void *, void *, void *, void *, void *, void*, void*);

/* .Call calls */
extern SEXP glm_deterministic(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_mcmc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_mcmcbas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_sampleworep(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP deterministic(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mcmc_new(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mcmcbas(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sampleworep_new(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);



static const R_CMethodDef CEntries[] = {
  {"gexpectations_vect", (DL_FUNC) &gexpectations_vect, 11},
  {"hypergeometric1F1",  (DL_FUNC) &hypergeometric1F1,  6},
  {"hypergeometric2F1",  (DL_FUNC) &hypergeometric2F1,  5},
  {"logHyperGauss2F1",   (DL_FUNC) &logHyperGauss2F1,   5},
  {"phi1",               (DL_FUNC) &phi1,               9},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"glm_deterministic", (DL_FUNC) &glm_deterministic, 11},
  {"glm_fit",           (DL_FUNC) &glm_fit,            7},
  {"glm_mcmc",          (DL_FUNC) &glm_mcmc,          16},
  {"glm_mcmcbas",       (DL_FUNC) &glm_mcmcbas,       16},
  {"glm_sampleworep",   (DL_FUNC) &glm_sampleworep,   15},
  {"mcmc_new",          (DL_FUNC) &mcmc_new,          20},
  {"deterministic",     (DL_FUNC) &deterministic,     11},
  {"mcmcbas",           (DL_FUNC) &mcmcbas,           17},
  {"sampleworep_new",   (DL_FUNC) &sampleworep_new,       15},
  {NULL, NULL, 0}
};

void R_init_BAS(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}



