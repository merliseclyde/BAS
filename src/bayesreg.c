#include "sampling.h"



void PrecomputeData(double *Xwork, double *Ywork, double *wts, double **pXtXwork, double **pXtYwork, double **pXtX, double **pXtY, double *yty, double *SSY, int p, int nobs) {
	char uplo[] = "U", trans[]="T";
	double one=1.0, zero=0.0;
	int inc=1, i, j,l;
	
	int p2 = p * p;
	*pXtXwork  = (double *) R_alloc(p2, sizeof(double));
	*pXtYwork = vecalloc(p);
	*pXtX  = (double *) R_alloc(p2, sizeof(double));
	*pXtY = vecalloc(p);

	for (j=0, l=0; j < p; j++) {
	  for	 (i=0; i < nobs; i++) {
	    Xwork[l] = Xwork[l]*wts[i];
	    l = l + 1; 
	  } 
	}
	
	for (i = 0; i< nobs; i++) {
	   Ywork[i] = Ywork[i]*wts[i];
	}
 
	//precompute XtX
	memset(*pXtX,0, p2 * sizeof(double));
	F77_NAME(dsyrk)(uplo, trans, &p, &nobs, &one, &Xwork[0], &nobs, &zero, *pXtX, &p); 

	double ybar = 0.0; 
	for (int i = 0; i< nobs; i++) {
		ybar += Ywork[i];
	}
	ybar = F77_NAME(ddot) (&nobs, &Ywork[0], &inc, &wts[0], &inc)/
	       F77_NAME(ddot) (&nobs, &wts[0], &inc, &wts[0], &inc);

	*yty = F77_NAME(ddot)(&nobs, &Ywork[0], &inc, &Ywork[0], &inc);	//	ybar = ybar/ (double) nobs;
	*SSY = *yty - (double) nobs* ybar *ybar;
	F77_NAME(dgemv)(trans, &nobs, &p, &one, &Xwork[0], &nobs, &Ywork[0], &inc, &zero, *pXtY,&inc);
}


/* Version of gexpectations that can be called directly from R */
void gexpectations_vect(int *nmodels, int *p, int *pmodel, int *nobs, double *R2, double *alpha, int *method,
                        double *RSquareFull, double *SSY, double *logmarg, double *shrinkage) {

  int i; 
     for (i=0; i < *nmodels; i++) {
	 gexpectations(*p, pmodel[i],  *nobs,  R2[i], *alpha, *method, 
		        *RSquareFull,  *SSY,  &logmarg[i],  &shrinkage[i]);

	     }
}



void gexpectations(int p, int pmodel, int nobs, double R2, double alpha, int method,  double RSquareFull, double SSY, double *logmarg, double *shrinkage) {  
    
    *shrinkage = 1.0;

    switch (method) { 
     case 0: 
      *logmarg = logBF_gprior(R2, nobs, pmodel,alpha); 
      if (pmodel > 1) *shrinkage = alpha/(1.0 + alpha);

      break;
     case 1:  
      *logmarg = logBF_hyperGprior(R2, nobs, pmodel,alpha); 
      *shrinkage = shrinkage_hyperg(R2, nobs, pmodel, alpha);
      break;
     case 2: 
      *logmarg = logBF_EB(R2, nobs, pmodel,alpha); 
      *shrinkage = shrinkage_EB_local(R2, nobs, pmodel,alpha);
      break;
     case 3: 
     *logmarg = BIC(R2, nobs, pmodel,SSY);
     *shrinkage = 1.0;
      break;
     case 4: 
       *logmarg = LogBF_ZS_null(R2, nobs, pmodel);
       *shrinkage = E_ZS_approx_null(R2,nobs,pmodel-1);
      break; 
     case 5: 
       *logmarg = LogBF_ZS_full(RSquareFull, R2, nobs, p, pmodel);
       *shrinkage = 1.0;
      break; 
      case 6: 
        *logmarg =  logBF_hyperGprior_laplace(R2, nobs, pmodel,alpha); 
	*shrinkage = shrinkage_laplace(R2, nobs, pmodel, alpha); 
      break; 
      case 7: 
	*logmarg =  AIC(R2,  nobs, pmodel,SSY);
	*shrinkage = 1.0;
      break; 
      case 8: 
	  *logmarg =  LogBF_Hg_null(R2,  nobs, pmodel, alpha, 1);
 	  if (pmodel > 1) {
	      *shrinkage = LogBF_Hg_null(R2,  nobs, pmodel+2, alpha, 2);
	      *shrinkage = exp(*shrinkage - *logmarg); }
      break; 
  default: 
      Rprintf("Error: Method must be one of g-prior, hyper-g, laplace (hyper-g), AIC, BIC, ZS-null, or ZS-full\n");
      break;
    }
}


double maxeffect(double *beta, double *se, int p) {
  int j;
  double t, maxeffect;
  
  maxeffect = 0;
  if (p > 1) {
    for (j = 2; j < p; j++) {
      t = fabs(beta[j]/se[j]);	
      maxeffect = fmax2(maxeffect, t);   
    }
  }
  return(maxeffect);
}	

// [[register]]
void cholreg(double *XtY, double *XtX, double *coefficients, double *se, double *mse, int p, int n)
{
	/* On entry *coefficients equals X'Y, which is over written with the OLS estimates */
	/* On entry MSE = Y'Y */

  double   ete, one, zero;
	int  job, l, i, j, info, inc;
	zero = 0.0;
	one = 1.0;
	inc = 1;
	job = 01;
	info = 1;


	/* LINPACK	
	F77_NAME(dpofa)(&XtX[0],&p,&p, &info);
	F77_NAME(dposl)(&XtX[0],&p,&p,&coefficients[0]);
	F77_NAME(dpodi)(&XtX[0],&p,&p, &det, &job);
	*/

	// LAPACK Equivalent
	//	F77_NAME(dpstrf)("U", &p, &XtX[0],&p, &piv[0], &rank, &tol, &work, &info);	

	F77_NAME(dpotrf)("U",&p, &XtX[0],&p, &info);
	F77_NAME(dpotrs)("U", &p, &inc, &XtX[0],&p,&coefficients[0],&p, &info);
	F77_NAME(dpotri)("U", &p, &XtX[0],&p, &job);
	ete = F77_NAME(ddot)(&p, &coefficients[0], &inc, &XtY[0], &inc);

	if ( n <= p) {
	  *mse = 0.0;
	}
	else {
	  *mse = (*mse - ete)/((double) (n - p));
	  }
	
	for (j=0, l=0; j < p; j++)  {
		for (i=0; i <  p; i++) {
			if (i == j)  {
				se[j] = sqrt(XtX[l] * *mse);
			}
			l += 1;    
		}
	}
} 


// [[register]]
double logBF_gprior( double Rsquare, int n,  int p, double g)
  {  double logmargy;
  /* log Marginal under reference prior for phi, intercept and
     g prior for regression coefficients 
     log marginal for null model is 0 */
    logmargy =  .5*(log(1.0 + g) * ((double) (n - p))  - log(1.0 + g * (1.0- Rsquare)) * ((double)(n - 1)));
  if (p == 1) logmargy = 0.0;
  return(logmargy);
  }

// [[register]]
double BIC(double Rsquare, int n,  int p, double SSY)
  {  double logmargy, sigma2, dn, dp;
  dn = (double) n;
  dp = (double) p;
  sigma2 = SSY*(1.0 - Rsquare);
  logmargy =  -.5*(dn*log(sigma2) +  dp*log(dn));
  return(logmargy);
  }

// [[register]]
double AIC(double Rsquare, int n,  int p, double SSY)
  {  double logmargy, sigma2, dn, dp;
  dn = (double) n;
  dp = (double) p;
  sigma2 = SSY*(1.0 - Rsquare);
  logmargy =  -.5*(dn*log(sigma2) +  dp*2.0);
  return(logmargy);
  }

// [[register]]
double shrinkage_EB_local(double R2, int n, int p, double alpha)
 {  
   /* R2 = Y'PxY/Y'Y  
      n = sample size
      p = number of rank of X (including intercept)
      alpha = prior hyperparameter
   */

   double  ghat, dn, dp, shrinkage;
   
    dn = (double) n - 1.0;
    dp = (double) p - 1.0;
    if (dp > 0.0) {
      ghat =  (((dn-dp)/dp)*R2/(1.0 - R2)) - 1.0;
      if (ghat < 0.0) ghat = 0.0;
      shrinkage = ghat/(1.0 + ghat);
      if (dp >= dn) shrinkage = 0.0;
    }
    else shrinkage = 1.0;
    
  return(shrinkage);
}

extern double hyp2f1(double a, double b, double c, double x);

double shrinkage_hyperg(double R2, int n, int p, double alpha)
{ double bnum, bden, a, c, Eg, Egplus1, z,s;
 a = ((double) n - 1.0 ) / 2.0;
 bnum = 2.0;
 bden = 1.0;
 c = (double) p  - 1.0 + alpha;
 if ( p == 1 || 2.0*a - c < 0.0) s = 1.0;
 else{ 	
   z = R2;
   Eg = hyp2f1(a, bnum, (c + 2.0)/2.0, z);
   Egplus1 =  hyp2f1(a, bden, c/2.0, z);
   s = 2.0*(Eg/Egplus1)/c;
   if (! R_FINITE(s)) s = shrinkage_laplace(R2, n, p, alpha);
 }
 return(s); 
}

// [[register]]
double logBF_hyperGprior(double R2, int n, int p, double alpha)
  {  double logmargy,  a, b, c, z1, hf1;
  /* log Marginal under reference prior for phi & intercept
     hyperg prior for regression coefficients; alpha > 2 
     log marginal for null model is 0 */

  a = (double)  (n - 1) /2.0;
  b = 1.0;
  c = ((double) p - 1.0  + alpha )/2.0;
  z1 = R2;
  logmargy = 0.0;
  if (a - c > 0) {
    hf1 = hyp2f1(a, b, c, z1);
    if (p == 1 || p >= n) logmargy = 0.0;
    else logmargy = log(hf1) 
	     - log( (double) p - 1.0 + alpha - 2.0) + log(2.0) 
	     + log(alpha/2.0 - 1.0);
    if (! R_FINITE(logmargy))
	logmargy = logBF_hyperGprior_laplace(R2, n, p, alpha);
  }
    return(logmargy);
  }

// [[register]]
double logBF_hyperGprior_laplace(double R2, int n, int p, double alpha)
 {  
   /* R2 = usual coefficient of determination
      n = sample size
      p = number of rank of X (including intercept)
      alpha = prior hyperparameter
      n and p are adjusted by subtrating off one 
      because of the flat prior on the intercept
   */

   double lognc, ghat, dn, dp, logmarg,sigmahat;
   
    dn = (double) n - 1.0;
    dp = (double) p - 1.0;
/*  Laplace approximation in terms of exp(tau) = g  */
/*  Agrees with Mathematica but not paper  */
  if (p == 1 || dn <= dp) logmarg = 0.0;
  else {
    ghat = (-4.+ alpha + dp + (2. - dn)*R2 - 
	    sqrt(-8.*(-2. + alpha + dp)*(-1.0 + R2) + (-4. + alpha + dp + (2.-dn)* R2)*(-4. + alpha + dp + (2.-dn)* R2)))/(2.*(-2. + alpha + dp)*(-1. + R2)); 

    if (ghat <= 0.0)  { Rprintf("ERROR: In Laplace approximation to  logmarg,  ghat =  %f  R2 = %f p = %d  n = %d\n", ghat, R2, p,n);}
  
   
    /*  Substitute ghat (unconstrained) into sigma, simplify, then plug in ghat
	Seems to perform better */

        
    sigmahat =1.0/(-ghat*(dn - alpha - dp)/(2.*(1. + ghat)*(1.+ghat)) +
                    dn*(ghat*(1. - R2))/(2.*(1.+ghat*(1.-R2))*(1.+ghat*(1.-R2)))); 

    if (sigmahat <= 0 ) Rprintf("ERROR in LAPLACE APPROXIMATION to logmarg sigmhat = %f, ghat =  %f  R2 = %f p = %d  n = %d\n", sigmahat, ghat, R2, p,n); 
    lognc = log(alpha/2.0 - 1.0);
    logmarg = lognc + 
              .5*( log(2.0*PI) 
                     - (dp + alpha)*log(1.0 + ghat)
	             -  dn*log(1.0-(ghat/(1.0 + ghat))*R2)
	             + log(sigmahat)) + log(ghat);
  }
  return(logmarg);
}

// [[register]]
double shrinkage_laplace(double R2, int n, int p, double alpha)
 {  
   /* R2 = usual R2
      n = sample size
      p = number of rank of X (including intercept)
      alpha = prior hyperparameter
   */

     double  ghat, dn, dp, lognum, logden, shrinkage,sigmahat, lognc;
   
   /* numerator Laplace approximation E[g/(1 + g) | Y, model] */

    dn = (double) (n - 1);
    dp = (double) (p - 1);
   if (p == 1) shrinkage = 1.0;
   else {
     if ( p >= n) shrinkage = 2.0/(2.0 + alpha);
     else {
       lognc = log(alpha/2.0 - 1.0);
       ghat = (-6.+ alpha + dp + (4. - dn)*R2 - 
	    sqrt(-16.*(-2. + alpha + dp)*(-1.0 + R2) + (-6. + alpha + dp + (4.-dn)* R2)*(-6. + alpha + dp + (4.-dn)* R2)))/(2.*(-2. + alpha + dp)*(-1. + R2)); 

       if (ghat <=0.0) { Rprintf("ERROR: In Laplace approximation to  E[g/(1 + g)] ghat = %f %f %d %d\n", ghat, R2, p,n ); 

       } 
       sigmahat =2.0/(ghat*(-dn + 2.+ alpha + dp)/((1. + ghat)*(1.+ghat)) +
		  dn*(ghat*(1. - R2))/((1.+ghat*(1.-R2))*(1.+ghat*(1.-R2)))); 

       if (sigmahat <= 0 ) Rprintf("ERROR in LAPLACE APPROXIMATION to E[g/(1 + g)] sigmahat = %f %f %f %d %d\n", sigmahat, ghat, R2, p,n); 

       lognum= .5*( log(2.0*PI) + 2.0*log(ghat)
		- (dp + alpha + 2.0 - dn)*log(1.0 + ghat)
		- dn*log(1.0 + ghat*(1. -R2))
		+ log(sigmahat))  +  lognc + log(ghat);
       logden = logBF_hyperGprior_laplace( R2,  n, p, alpha);
   // lognc is included here so that it cancels wth lognc for denominator
       shrinkage = exp(lognum - logden);
   //  Rprintf("%f %f %f  %f %f %f \n", ghat, sigmahat, normalpart,
   //  lognum, logden, shrinkage);
     }}
  return(shrinkage);
 }

// [[register]]
double log_laplace_2F1(double a, double b, double c, double z)
 {  
   
   double ghat1, ghat2, ghat,logint,sigmahat, D, A, B, C,root1, root2;

   // 2F1(a,b,c, Z)      =        Gamma(c)
   //                        --------------------- * integral
   //			     Gamma(b) Gamma(c - b)

   // integral = \int_0^\infty  g^{b -1} (1 + g)^{c - b - 1}(1 + (1 - z)g)^{-a}

   logint = 0.0;

   if (c >= b && b > 0) {
     logint = lgammafn(c) - lgammafn(b) - lgammafn(c - b);
   }
   
   if ( z == 1.0) {
     if (c > b + a) 
     logint =  lgammafn(c) + lgammafn(c - a - b) - lgammafn(c - b)
       - lgammafn(c - a);
     else logint = log(0.0);
   }
   else{	
/*  Laplace approximation in terms of exp(tau) = g  */
//
//  Integrate[g^(b-1) (1 + g)^(a - c) (1 + (1 - z) g)^(-a) dg 
//  Integrate[e^tau b (1 + e^g)^(a - c) *( 1 + (1 -z)e^g)^(-a) dtau  
//  Jacobian e^tau 
   
   A = (b - c)*(1. - z);
   B = 2.*b - c + (a-b)*z;
   C = b;
   D = B*B - 4.*A*C;

   if (D < 0 )    Rprintf("ERROR: In Laplace approximation to 2F1");

   root1 = (-B - sqrt(D))/(2.*A);
   root2 = (-B + sqrt(D))/(2.*A);
   
   ghat1 = (2.*b)/(c + b*(z - 2.0) - a*z + sqrt(c*c + 4*a*b*z -  2.0*(a + b)*c *z + (a -b )*(a-b)*z*z));
 
   ghat2 =  -(c +b*(-2. + z) - a*z + sqrt(c*c + 4*a*b*z -  2.0*(a + b)*c *z + (a -b )*(a-b)*z*z))/(2.*(b - c)*(z - 1.));
   
   // Rprintf("+root = %lf -root = %lf ghat1  = %lf hgat2 = %lf\n", root2, root1, ghat1, ghat2);
   ghat = ghat1;
   if ( ghat1 < 0) Rprintf("ERROR: In Laplace approximation to 2F1");

   sigmahat =1.0/((-a + c)*((ghat/(1 + ghat))*(1 - ghat/(1 + ghat))) +
		    a*((ghat*(1.-z))/(1.+ghat*(1.-z))*
		       (1. - (ghat*(1.-z))/(1.+ghat*(1.-z)))));

    if (sigmahat <= 0 ) Rprintf("ERROR in LAPLACE APPROXIMATION to in 2F1 sigmhat = %f, ghat =  %f  z = %f \n", sigmahat, ghat, z); 
    logint = logint
                   + .5*( log(2.0*PI) +  log(sigmahat)) +
		      b*log(ghat) + (a - c)*log(1.0 + ghat)
                   -  a*log(1.0 + ghat*(1. - z)); 
   }
  return(logint);
 }	

// [[register]]
void logHyperGauss2F1(double *a, double  *b, double *c, double *z, double *integral){
  // version that can be called from R using .Call 
  *integral = log_laplace_2F1(*a, *b, *c, *z);
}

// [[register]]
double logBF_EB(double R2, int n, int p, double alpha)
 {
   double  ghat, dn, dp, logmarg;
    dn = (double) n -1.0;
    dp = (double) p -1.0;
    ghat =  (((dn-dp)/dp)*R2/(1.0 - R2)) - 1.0;
    if (ghat < 0.0) ghat = 0.0;
    logmarg = .5*( - dn*log(1.0-(ghat/(1.0 + ghat))*R2)
		   - dp*log(1.0 + ghat));
    if (p == 1 || p >= n) logmarg = 0.0;
    return(logmarg);
}

double CalculateRSquareFull(double *XtY, double *XtX, double *XtXwork, double *XtYwork, 
							SEXP Rcoef_m, SEXP Rse_m, int p, int nobs, double yty, double SSY) {
	double RSquareFull;
	if (nobs <= p) {
		RSquareFull = 1.0;
	} else {
		PROTECT(Rcoef_m = NEW_NUMERIC(p));
		PROTECT(Rse_m = NEW_NUMERIC(p));
		double *coefficients = REAL(Rcoef_m);  
		double *se_m = REAL(Rse_m);
		memcpy(coefficients, XtY,  p*sizeof(double));
		memcpy(XtXwork, XtX, p * p *sizeof(double));
		memcpy(XtYwork, XtY,  p*sizeof(double));

		double mse_m = yty; 
		cholreg(XtYwork, XtXwork, coefficients, se_m, &mse_m, p, nobs);  

		RSquareFull =  1.0 - (mse_m * (double) ( nobs - p))/SSY;
		UNPROTECT(2);
	}
	return RSquareFull;
}
