#' Empirical Bayes Prior Distribution for Coefficients in BMA Model
#' 
#' Creates an object representing the EB prior for BAS GLM.
#' 
#' Creates a structure used for \code{\link{bas.glm}}.
#' 
#' @aliases EB EB.local
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{CCH}} and \code{\link{bas.glm}}
#' @examples
#' EB.local()
#' 
#' 
#' @rdname EB.local
#' @family beta priors
#' @export

EB.local = function() {
    structure(list(family="EB-local", class="EB", hyper.parameters=list(local=TRUE)),
    class="prior")
}



#' Generalized g-Prior Distribution for Coefficients in BMA Models
#' 
#' Creates an object representing the CCH mixture of g-priors on coefficients
#' for BAS .
#' 
#' Creates a structure used for \code{\link{bas.glm}}.
#' 
#' @param alpha a scalar > 0, recommended alpha=.5 (betaprime) or 1 for CCH.
#' The hyper.g(alpha) is equivalent to CCH(alpha -2, 2, 0). Liang et al
#' recommended values in the range 2 < alpha_h <= 4
#' @param beta a scalar > 0.  The value is not updated by the data; beta should
#' be a function of n for consistency under the null model.  The hyper-g
#' corresonds to b = 2
#' @param s a scalar, recommended s=0
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise A Clyde
#' @seealso \code{\link{IC.prior}}, \code{\link{bic.prior}},
#' \code{\link{bas.glm}}
#' @examples
#' CCH(alpha=.5, beta=100, s=0) 
#' 
#' @rdname CCH
#' @family beta priors
#' @export
#' 
#' 
CCH = function(alpha, beta, s=0) {
#    if (beta == 2 & alpha == 2 & s == 0)   {
#        structure(list(family="Truncated-Gamma", class="TCCH", hyper.parameters=NULL),
#                       class="prior")}
#    else {      
        structure(list(family="CCH", class="TCCH",
                       hyper.parameters=list(alpha=alpha, beta=beta, s=s)),
                  class="prior")
    #}
}



#' Generalized tCCH g-Prior Distribution for Coefficients in BMA Models
#' 
#' Creates an object representing the tCCH mixture of g-priors on coefficients
#' for BAS.
#' 
#' Creates a structure used for \code{\link{bas.glm}}.
#' 
#' @param alpha a scalar > 0, recommended alpha=.5 (betaprime) or 1.
#' @param beta a scalar > 0.  The value is not updated by the data; beta should
#' be a function of n for consistency under the null model.
#' @param s a scalar, recommended s=0 a priori
#' @param r r arbitrary; in the hyper-g-n prior sets r = (alpha + 2)
#' @param v 0 < v
#' @param theta theta > 1
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{CCH}}, \code{\link{robust}}, \code{\link{hyper.g}},
#' \code{\link{hyper.g.n}}\code{\link{bas.glm}}
#' @examples
#' n = 500;
#'  tCCH(alpha=1, beta=2, s=0, r=1.5, v = 1, theta=1/n)
#' @rdname tCCH
#' @family beta priors
#' @export
tCCH = function(alpha=1, beta=2, s=0, r=3/2, v=1, theta=1) {
#    if (beta == 2 & alpha == 2 & s == 0)   {
#        structure(list(family="Truncated-Gamma", class="TCCH", hyper.parameters=NULL),
#                       class="prior")}
#    else {      
        structure(list(family="tCCH", class="TCCH",
                       hyper.parameters=list(alpha=alpha, beta=beta, s=s,
                           r=r, v=v, theta=theta)),
                  class="prior")
    #}
    }



#' Intrinsic Prior Distribution for Coefficients in BMA Models
#' 
#' Creates an object representing the intrinsic prior on g, a special case of
#' the tCCH mixture of g-priors on coefficients for BAS.
#' 
#' Creates a structure used for \code{\link{bas.glm}}.
#' 
#' @param n the sample size; if NULL, the value derived from the data in the
#' call to `bas.glm` will be used.
#' @return returns an object of class "prior", with the family "intrinsic" of
#' class "TCCH" and hyperparameters alpha = 1, beta = 1, s = 0, r = 1, n = n
#' for the tCCH prior where theta in the tCCH prior is determined by the model
#' size and sample size.
#' @author Merlise A Clyde
#' @seealso \code{\link{tCCH}}, \code{\link{robust}}, \code{\link{hyper.g}},
#' \code{\link{hyper.g.n}}\code{\link{bas.glm}}
#' @references Womack, A., Novelo,L.L., Casella, G. (2014). "Inference From
#' Intrinsic Bayes' Procedures Under Model Selection and Uncertainty". Journal
#' of the American Statistical Association.  109:1040-1053.
#' \url{http://amstat.tandfonline.com/doi/abs/10.1080/01621459.2014.880348}
#' 
#' @examples
#' n = 500;
#'  tCCH(alpha=1, beta=2, s=0, r=1.5, v = 1, theta=1/n)
#' 
#' 
#' @rdname intrinsic
#' @family beta priors
#' @export

intrinsic = function(n=NULL) {
#    if (beta == 2 & alpha == 2 & s == 0)   {
#        structure(list(family="Truncated-Gamma", class="TCCH", hyper.parameters=NULL),
#                       class="prior")}
#    else {      
        structure(list(family="intrinsic", class="TCCH",
                       hyper.parameters=list(alpha=1.0, beta=1.0, s=0.0, r=1.0, n=n)),
                  class="prior")
    #}
    }



#' Generalized hyper-g/n Prior Distribution for g for mixtures of g-priors on
#' Coefficients in BMA Models
#' 
#' Creates an object representing the hyper-g/n mixture of g-priors on
#' coefficients for BAS. This is a special case of the tCCH prior
#' 
#' Creates a structure used for \code{\link{bas.glm}}.  This is a special case
#' of the \code{\link{tCCH}}, where \code{hyper.g.n(alpha=3, n)} is equivalent
#' to \code{ tCCH(alpha=1, beta=2, s=0, r=1.5, v = 1, theta=1/n) }
#' 
#' @param alpha a scalar > 0, recommended 2 < alpha <= 3
#' @param n The sample size; if NULL, the value derived from the data in the
#' call to `bas.glm` will be used.
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{tCCH}}, \code{\link{robust}}, \code{\link{hyper.g}},
#' \code{\link{CCH}}\code{\link{bas.glm}}
#' @examples
#' n = 500
#' hyper.g.n(alpha = 3, n=n)
#'
#' @rdname hyper.g.n  
#' @family beta priors
#' @export
hyper.g.n = function(alpha=3, n=NULL) {
#    if (beta == 2 & alpha == 2 & s == 0)   {
#        structure(list(family="Truncated-Gamma", class="TCCH", hyper.parameters=NULL),
#                       class="prior")}
#    else {      
        structure(list(family="hyper-g/n", class="TCCH",
                       hyper.parameters=list(alpha=alpha-2, beta=2, s=0,
                           r=alpha/2, v=1, theta=1/n), n=n),
                  class="prior")
    #}
}



#' Jeffreys Prior Distribution for $g$ for Mixtures of g-Priors for
#' Coefficients in BMA Models
#' 
#' Creates an object representing the Jeffrey's Prior on g mixture of g-priors
#' on coefficients for BAS. This is equivalent to a limiting version of the
#' CCH(a, 2, 0) with a = 0 or they hyper-g(a = 2) and is an improper prior.  As
#' $g$ does not appear in the Null Model, Bayes Factors and model probabilities
#' are not well-defined because of arbitrary normalizing constants, and for
#' this reason the null model is excluded and the same c onstants are used
#' across other models.
#' 
#' Creates a structure used for \code{\link{bas.glm}}.
#' 
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{CCH}} \code{\link{bas.glm}}
#' @examples
#' Jeffreys()
#' 
#' 
#' @rdname Jeffreys
#' @family beta priors
#' @export
Jeffreys = function() {
        structure(list(family="Jeffreys", class="TCCH",
                       hyper.parameters=list(alpha=0, beta=2, s=0.0)),
                  class="prior")
    }




#' Hyper-g-Prior Distribution for Coefficients in BMA Models
#' 
#' Creates an object representing the hyper-g mixture of g-priors on
#' coefficients for BAS.
#' 
#' Creates a structure used for \code{\link{bas.glm}}.
#' 
#' @param alpha a scalar > 0. The hyper.g(alpha) is equivalent to CCH(alpha -2,
#' 2, 0). Liang et al recommended values in the range 2 < alpha_h <= 3
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{CCH}} \code{\link{bas.glm}}
#' @examples
#' hyper.g(alpha=.5) 
#' 
#' 
#' @rdname hyper.g
#' @family beta priors
#' @export
hyper.g = function(alpha=3.0) {
    if (alpha <= 2 )   {
        return("alpha must be greater than 2 in hyper.g prior")
    }
    
    else {      
        structure(list(family="CCH", class="TCCH",
                       hyper.parameters=list(alpha=alpha-2.0, beta=2, s=0.0)),
                  class="prior")}
}



#' Generalized g-Prior Distribution for Coefficients in BMA Models
#' 
#' Creates an object representing the Truncated Gamma (tCCH) mixture of
#' g-priors on coefficients for BAS, where u = 1/(1+g) has a Gamma distribution
#' supported on (0, 1].
#' 
#' Creates a structure used for \code{\link{bas.glm}}.
#' 
#' @param alpha a scalar > 0, recommended alpha=.5 (betaprime) or 1.  alpha=2
#' corresponds to the uniform prior on the shrinkage factor.
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{CCH}} \code{\link{bas.glm}}
#' @examples
#' 
#' TG(alpha=2)
#' CCH(alpha=2, beta=100, s=0)
#' 
#' @family beta priors
#' @export
TG = function(alpha=2) {
    structure(list(family="TG", class="TCCH",
                       hyper.parameters=list(alpha=alpha, beta=2.0, s=0.0)),
                  class="prior")
}




#' Beta-Prime Prior Distribution for Coefficients in BMA Model
#' 
#' Creates an object representing the Beta-Prime prior that is mixture of
#' g-priors on coefficients for BAS.
#' 
#' Creates a structure used for \code{\link{bas.glm}}.
#' 
#' @param n the sample size; if NULL, the value derived from the data in the
#' call to `bas.glm` will be used.
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{CCH}}
#' @examples
#' beta.prime(n=100)
#' 
#' @rdname beta.prime
#' @family beta priors
#' @export

beta.prime = function(n=NULL) {
    structure(list(family="betaprime", class="TCCH", 
                   hyper.parameters=list(n=n, alpha=.5)),
              class="prior")
}



#' Robust-Prior Distribution for Coefficients in BMA Model
#' 
#' Creates an object representing the robust prior of Bayarri et al (2012) that
#' is mixture of g-priors on coefficients for BAS.
#' 
#' Creates a prior structure used for \code{\link{bas.glm}}.
#' 
#' @param n the sample size.
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{CCH}} and\code{\link{bas.glm}}
#' @examples
#' robust(100)
#' 
#' 
#' @rdname robust
#' @family beta priors
#' @export

robust = function(n=NULL) {
    structure(list(family="robust", class="TCCH",  
                   hyper.parameters = list(n=n)),
                   class="prior")
}

#' @family prior functions
#' @export

bic.prior = function(n=NULL) {
   if (is.null(n)) penalty = NULL
   else penalty = log(n)
   
  structure(list(family="BIC", class="IC", 
                   hyper.parameters = list(penalty=penalty, n=as.numeric(n)),
                   hyper=penalty),
                   class="prior")
}

#' @family beta priors
#' @export

aic.prior = function() {
    structure(list(family="AIC", class="IC", hyper.parameters = list(penalty=2),
                   hyper=2),
                   class="prior")
}

    


#' Information Criterion Families of Prior Distribution for Coefficients in BMA
#' Models
#' 
#' Creates an object representing the prior distribution on coefficients for
#' BAS.
#' 
#' The log marginal likelihood is approximated as -2*(deviance +
#' penalty*dimension).  Allows alternatives to AIC (penalty = 2) and BIC
#' (penalty = log(n)).  For BIC, the argument may be missing, in which case the
#' sample size is determined from the call to `bas.glm` and used to dertermine
#' the penalty.
#' 
#' @aliases IC.prior aic.prior AIC.prior bic.prior BIC.prior
#' @param penalty a scalar used in the penalized loglikelihood of the form
#' penalty*dimension
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{g.prior}}
#' @examples
#' IC.prior(2)
#'           aic.prior()
#'           bic.prior(100)
#'           
#' @family beta priors 
#' @export

IC.prior = function(penalty) {
    structure(list(family="IC", class="IC", hyper=as.numeric(penalty),
                   hyper.parameters = list(penalty=penalty)),
              class="prior")
}    


#' Families of G-Prior Distribution for Coefficients in BMA Models
#' 
#' Creates an object representing the g-prior distribution on coefficients for
#' BAS.
#' 
#' Creates a structure used for BAS.
#' 
#' @param g a scalar used in the covariance of Zellner's g-prior, Cov(beta) =
#' sigma^2 g (X'X)^-1
#' 
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{IC.prior}}
#' @examples
#' g.prior(100)
#' 
#' @family beta priors
#' @export

g.prior = function(g) {
    structure(list(family="g.prior", g = as.numeric(g), class="g-prior",
                   hyper=as.numeric(g),
                   hyper.parameters= list(g=g)),
              class="prior")
}



#' Test based Bayes Factors for BMA Models
#' 
#' Creates an object representing the prior distribution on coefficients for
#' BAS that coreesponds to the test-based Bayes Factors.
#' 
#' Creates a prior object structure used for BAS in `bas.glm`.
#' 
#' @param g a scalar used in the covariance of Zellner's g-prior, Cov(beta) =
#' sigma^2 g (X'X)^-
#' 
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{g.prior}}, \code{\link{bas.glm}}
#' @examples
#' 
#' testBF.prior(100)
#' library(MASS)
#' data(Pima.tr)
#' 
#' # use g = n 
#' bas.glm(type ~ ., data=Pima.tr, family=binomial(), 
#'         betaprior=testBF.prior(nrow(Pima.tr)),
#'         modelprior=uniform(), method="BAS")
#' @family beta priors 
#' @export

testBF.prior = function(g) {
  structure(list(family="testBF.prior", g = as.numeric(g), class="g-prior",
                 hyper=as.numeric(g),
                 hyper.parameters= list(g=as.numeric(g), loglik_null=NULL)),
            class="prior")
}
