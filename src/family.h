/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2005   The R Development Core Team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

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
  double (*info_matrix) (double *y, double *eta, int n);
} glmstptr;

glmstptr * make_glmfamily_structure(SEXP family);

// Binomial family
double binomial_loglik(double *Y, double *mu, int n);
void logit_link(double *mu, double *eta, int n);
void logit_variance(double *mu, double *var, int n);
void logit_linkinv(double *eta, double *mu, int n);
void logit_mu_eta(double *eta, double *mu, int n);
void logit_info(double *y, double *mu, double *Imu, int n);
void binomial_dev_resids(double *y, double *mu, double *wt, double *res, int n);
double binomial_dispersion(double *resid,  double *weights, int n, int rank); 
void binomial_initialize(double *Y, double *mu,  double *weights, int n);
// Poisson
double poisson_dispersion(double *resid,  double *weights, int n, int rank) ;
void log_link(double *mu, double *eta, int n);
void poisson_variance(double *mu, double *var, int n);
void log_linkinv(double *eta, double *mu, int n);
void log_mu_eta(double *eta, double *mu, int n);
void poisson_log_info(double *y, double *mu, double *Imu, int n);
void poisson_dev_resids(double *y, double *mu, double *wt, double *res, int n);
void poisson_initialize(double *Y, double *mu,  double *weights, int n);


// Normal
double Gaussian_dispersion(double *resid,  double *weights, int n, int rank) ;



// General Functions
double deviance(double *res, int n);
void chol2se(double *qr, double *se, double *cov, double *covwork, int p, int n);
void QR2cov(double *qr, double *cov, double *covwork, int p, int n);
#endif

// Marginal likelihood, shrinkage and posterior expctation functions 

