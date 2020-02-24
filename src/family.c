/* Code for implementing GLM families in C   */

# include "bas.h"

static const double THRESH = 30.;
static const double MTHRESH = -30.;
static const double INVEPS = 1/DOUBLE_EPS;

/**
 * Evaluate x/(1 - x). An inline function is used so that x is
 * evaluated once only.
 *
 * @param x input in the range (0, 1)
 *
 * @return x/(1 - x)
 */
static R_INLINE double x_d_omx(double x) {
    if (x < 0 || x > 1)
	error(_("Value %d out of range (0, 1)"), x);
    return x/(1 - x);
}

/**
 * Evaluate x/(1 + x). An inline function is used so that x is
 * evaluated once only. [but inlining is optional!]
 *
 * @param x input
 *
 * @return x/(1 + x)
 */
static R_INLINE double x_d_opx(double x) {return x/(1 + x);}


double poisson_loglik(double *Y, double*mu, double *wts, int n) {
  int i;
  double ll = 0.0;

  for (i = 0; i < n; i++) {
    ll += wts[i]*dpois(Y[i],mu[i],1);
  }
  return(ll);
}


void poisson_variance(double *mu, double *var, int n) {

  int i;

  for (i = 0; i<n; i++) {
    var[i] = mu[i];
  }
}


void poisson_log_info(double *y, double *mu, double *weights, double *var, int n) {

  int i;
//  CHECK IF USE OF WEIGHTS IS CORRECT
  for (i = 0; i<n; i++) {
    var[i] = weights[i]*mu[i];
  }
}


void log_link(double *rmu, double *rans, int n) {
    int i;

    for (i = 0; i < n; i++)
	rans[i] = log(rmu[i]);
}


void log_linkinv(double *reta, double *rans, int n) {

    int i;

    for (i = 0; i < n; i++) {
      rans[i] = fmax2(exp(reta[i]), DOUBLE_EPS);
    }
}


void log_mu_eta(double *reta, double *rans, int n)
{
    int i;
    for (i = 0; i < n; i++) {
      rans[i] = fmax2(exp(reta[i]), DOUBLE_EPS);
    }
}


void poisson_dev_resids(double *ry, double *rmu, double *rwt, double *rans, int n)
{
  int i;
  double mui, yi, wti;

  	for (i = 0; i < n; i++) {
	    mui = rmu[i];
	    yi = ry[i];
	    wti = rwt[i];
	    rans[i] = mui * wti;
	    if (yi > 0) {
	      rans[i] = wti*(yi*log(yi/mui) - (yi - mui));
	    }
	    rans[i] *= 2.0;
	}
}


void poisson_initialize(double *Y, double *mu,  double *weights, int n) {
  int i;
  for (i = 0; i < n; i++) {
    if (Y[i] < 0.0) error("negative values not allowed for Poisson");
    mu[i] =  Y[i] + 0.1;
  }
}


double poisson_dispersion(double *resid,  double *weights, int n, int rank) {
  return(1.0);
}


/* Gamma */

double gamma_loglik(double *Y, double*mu, double *wts, int n) {
  int i;
  double ll = 0.0;
  double dev = 1.0, disp; 

  disp = dev/n;
  
  for (i = 0; i < n; i++) {
    ll += wts[i]*dgamma(Y[i],1/disp,1/(mu[i]*disp),1);
  }
  return(ll);
}


void gamma_variance(double *mu, double *var, int n) {
  
  int i;
  
  for (i = 0; i<n; i++) {
    var[i] = pow (mu[i],2.0);
  }
}


void gamma_dev_resids(double *ry, double *rmu, double *rwt, double *rans, int n)
{
  int i;
  double mui, yi, wti;
  
  for (i = 0; i < n; i++) {
    mui = rmu[i];
    yi = ry[i];
    wti = rwt[i];
    rans[i] = 2.0 *  wti * (yi - mui)/mui;
    if (yi > 0) {
      rans[i] += -2.0 * wti * log(yi/mui) ;
    }
  }
}


void gamma_initialize(double *Y, double *mu,  double *weights, int n) {
  int i;
  for (i = 0; i < n; i++) {
    if (Y[i] < 0.0) error("negative values not allowed for Gamma");
    mu[i] =  Y[i];
  }
}


double gamma_dispersion(double *resid,  double *weights, int n, int rank) {
  return(1.0);
}



/* Binomial */

double binomial_loglik(double *Y, double*mu, double *wts, int n) {
  int i;
  double ll = 0.0;

  for (i = 0; i < n; i++) {
    ll += wts[i]*dbinom(Y[i],1.0,mu[i],1);
  }
  return(ll);
}


void logit_variance(double *mu, double *var, int n) {

  int i;

  for (i = 0; i<n; i++) {
    var[i] = mu[i]*(1.0 - mu[i]);
  }
}


void logit_info(double *y, double *mu, double *weights, double *var, int n) {

  int i;

  for (i = 0; i<n; i++) {
    var[i] = mu[i]*(1.0 - mu[i])*weights[i];
  }
}


void logit_precision(double *mu, double *prec, int n) {

  int i;

  for (i = 0; i<n; i++) {
    prec[i] = 1/(mu[i]*(1.0 - mu[i]));
  }
}


void logit_link(double *rmu, double *rans, int n)
{
    int i;

    for (i = 0; i < n; i++)
	rans[i] = log(x_d_omx(rmu[i]));
}


void logit_linkinv(double *reta, double *rans, int n) {

    int i;

    for (i = 0; i < n; i++) {
	double etai = reta[i], tmp;
	tmp = (etai < MTHRESH) ? DOUBLE_EPS :
	    ((etai > THRESH) ? INVEPS : exp(etai));
	rans[i] = x_d_opx(tmp);
    }
}


void logit_mu_eta(double *reta, double *rans, int n)
{
    int i;
    for (i = 0; i < n; i++) {
		double etai = reta[i];
		double opexp = 1 + exp(etai);

		rans[i] = (etai > THRESH || etai < MTHRESH) ? DOUBLE_EPS :
			exp(etai)/(opexp * opexp);
    }
}

static R_INLINE
double y_log_y(double y, double mu)
{
    return (y) ? (y * log(y/mu)) : 0;
}


void binomial_dev_resids(double *ry, double *rmu, double *rwt, double *rans, int n)
{
  int i;
  double mui, yi;
  	for (i = 0; i < n; i++) {
	    mui = rmu[i];
	    yi = ry[i];
	    rans[i] = 2 * rwt[i] *
		(y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
	}
}


double binomial_dispersion(double *resid,  double *weights, int n, int rank) {
  return(1.0);
}


void binomial_initialize(double *Y, double *mu,  double *weights, int n) {
  int i;
  for (i = 0; i < n; i++) {
    if (weights[1] == 0) Y[i] = 0.0;
    mu[i] = (weights[i] * Y[i] + 0.5)/(weights[i] + 1.0) ;
  }
}


/* Gaussian */

double Gaussian_dispersion(double *resid,  double *weights, int n, int rank) {
  double dispersion = 0.0;
  int i, nwt=0;

  for (i = 0; i<n; i++) {
    if (weights[i] > 0) nwt += 1;
    dispersion += weights[i]*resid[i]*resid[i];
  }
   return(dispersion/(double) (nwt - rank));
}


/* generic functions */


double deviance(double *res, int n) {
  int i;
  double   dev = 0;

  for (i=0; i<n; i++) {
    dev += res[i];
  }
  return dev;
}


double quadform (double *bwork, double *R,  int p) {

  double Q = 0.0;
  int inc = 1;
  char uplo[] = "U", trans[]="T", diag[]="N";
  //  F77_NAME(dcopy)(&p, &b[0], &inc,  &bwork[0], &inc);
  F77_NAME(dtrsv)(uplo, trans, diag, &p, &R[0], &p, &bwork[0], &inc);
  Q = F77_NAME(dnrm2)(&p, &bwork[0], &inc);
  Q *=Q;
  return(Q);
}


void chol2se(double *qr, double *se, double *R, double *covwork, int p, int n) {

  int i, j, l;

  for (j=0, l=0; j < p; j++) {

    for (i = 0; i <p; i++, l++) {
      R[l] = 0;
      if (i < (j+1))   R[j*p+i] = qr[j*n + i];
    }
  }

  Lapack_chol2inv(R, p, covwork);

for (j=0; j < p; j++) {
  se[j] = sqrt(covwork[j*p + j]);
}

 return;
}


void QR2cov(double *qr,  double *R, double *covwork, int p,  int n) {

  int i, j, l;

  for (j=0, l=0; j < p; j++) {

    for (i = 0; i <p; i++, l++) {
      R[l] = 0;
      if (i < (j+1))   R[j*p+i] = qr[j*n + i];
    }
  }

  Lapack_chol2inv(R, p, covwork);
  return;
}



void  Lapack_chol2inv(double *A, int sz, double *ans)
{
  int  i, j;
  //	F77_NAME(dcopy)(&sz, &A[0], &inc,  &ans[0], &inc);
	for (j = 0; j < sz; j++) {
	    for (i = 0; i <= j; i++)
		ans[i + j * sz] = A[i + j * sz];
	}

	F77_CALL(dpotri)("Upper", &sz, &ans[0], &sz, &i);
	if (i != 0) {
	    if (i > 0)
		error(_("element (%d, %d) is zero, so the inverse cannot be computed"),
		      i, i);
	    error(_("argument %d of Lapack routine %s had invalid value"),
		  -i, "dpotri");
	}

	for (j = 0; j < sz; j++) {
	    for (i = j+1; i < sz; i++)
		ans[i + j * sz] = ans[j + i * sz];
	}
}






struct glmfamilystruc * make_glmfamily_structure(SEXP family) {

  glmstptr *glmfamily;

  //  Rprintf("Make family\n");

  glmfamily = (struct glmfamilystruc *) R_alloc(1, sizeof(struct glmfamilystruc));
  glmfamily->family = CHAR(STRING_ELT(getListElement(family, "family"),0));
	//	Rprintf("family %s\n", glmfamily->family);

	glmfamily->link = CHAR(STRING_ELT(getListElement(family, "link"),0));
	//	Rprintf("link %s\n", glmfamily->link);

	if  (strcmp(glmfamily->family, "binomial") == 0) {
		glmfamily->dev_resids = binomial_dev_resids;
		glmfamily->dispersion = binomial_dispersion;
		glmfamily->initialize = binomial_initialize;
		glmfamily->loglik = binomial_loglik;
		if (strcmp(glmfamily->link, "logit") != 0) {
		  warning("no other links implemented yet, using logit\n");
		}

		glmfamily->linkfun = logit_link;
		glmfamily->mu_eta = logit_mu_eta;
		glmfamily->variance = logit_variance;
		glmfamily->linkinv =  logit_linkinv;
		glmfamily->info_matrix =  logit_info;
	}
	else if  (strcmp(glmfamily->family, "poisson") == 0) {
	  glmfamily->dev_resids = poisson_dev_resids;
	  glmfamily->dispersion = poisson_dispersion;
	  glmfamily->initialize = poisson_initialize;
	  glmfamily->variance = poisson_variance;
	  glmfamily->loglik =  poisson_loglik;
	  if (strcmp(glmfamily->link, "log") != 0) {
	    warning("no other links implemented yet, using log\n");
	  }
	  glmfamily->linkfun = log_link;
	  glmfamily->mu_eta = log_mu_eta;
	  glmfamily->linkinv =  log_linkinv;
	  glmfamily->info_matrix =  poisson_log_info;
	}
	else if  (strcmp(glmfamily->family, "Gamma") == 0) {
	  glmfamily->dev_resids = gamma_dev_resids;
	  glmfamily->dispersion = gamma_dispersion;
	  glmfamily->initialize = gamma_initialize;
	  glmfamily->variance = gamma_variance;
	  glmfamily->loglik =  gamma_loglik;
	  if (strcmp(glmfamily->link, "log") != 0) {
	    warning("no other links implemented yet, using log\n");
	  }
	  glmfamily->linkfun = log_link;
	  glmfamily->mu_eta = log_mu_eta;
	  glmfamily->linkinv =  log_linkinv;
	  glmfamily->info_matrix =  poisson_log_info;
	}

	else {
	  error("only 'binomial() and 'poisson() and 'gamma() families supported now\n");
	  //	  stop(1);
	}
	return(glmfamily);
}
