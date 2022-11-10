#' Fitting Generalized Linear Models and Bayesian marginal likelihood
#' evaluation
#'
#' A version of glm.fit rewritten in C; also returns marginal likelihoods for
#' Bayesian model comparison
#'
#' C version of glm-fit.  For different prior choices returns, marginal
#' likelihood of model using a Laplace approximation.
#' @rdname bayesglm.fit
#' @param x design matrix
#' @param y response
#' @param weights optional vector of weights to be used in the fitting process.
#' Should be NULL or a numeric vector.
#' @param start starting value for coefficients in the linear predictor
#' @param etastart starting values for the linear predictor
#' @param mustart starting values for the vectors of means
#' @param offset a priori known component to be included in the linear
#' predictor
#' @param family a description of the error distribution and link function for
#' exponential family; currently only binomial(), poisson(), and Gamma() with canonical
#' links are implemented.
#' @param coefprior function specifying prior distribution on coefficients with
#' optional hyperparameters leading to marginal likelihood calculations;
#' options include \code{bic.prior()},\code{ aic.prior()}, and
#' \code{ic.prior()}
#' @param control a list of parameters that control convergence in the fitting
#' process.  See the documentation for \code{glm.control()}
#' @param intercept should an intercept be included in the null model?
#' @return \item{coefficients}{MLEs} \item{se}{Standard errors of coefficients
#' based on the sqrt of the diagonal of the inverse information matrix}
#' \item{mu}{fitted mean} \item{rank}{numeric rank of the fitted linear model}
#' \item{deviance}{minus twice the log likelihood evaluated at the MLEs}
#' \item{g}{value of g in g-priors} \item{shrinkage}{shrinkage factor for
#' coefficients in linear predictor} \item{RegSS}{quadratic form
#' beta'I(beta)beta used in shrinkage} \item{logmarglik}{the log marginal or
#' integrated log likelihood (up to a constant)}
#' @author Merlise Clyde translated the \code{\link{glm.fit}} from R base into
#' C using the .Call interface
#' @seealso \code{\link{bic.prior}}
#' @references \code{\link{glm}}
#' @keywords regression GLM
#' @examples
#' data(Pima.tr, package="MASS")
#' Y <- as.numeric(Pima.tr$type) - 1
#' X <- cbind(1, as.matrix(Pima.tr[,1:7]))
#' out <- bayesglm.fit(X, Y, family=binomial(),coefprior=bic.prior(n=length(Y)))
#' out$coef
#' out$se
#' # using built in function
#' glm(type ~ ., family=binomial(), data=Pima.tr)
#'
#' @export
#'
bayesglm.fit <-
  function(x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
             mustart = NULL, offset = rep(0, nobs), family = binomial(),
             coefprior = bic.prior(nobs),
             control = glm.control(), intercept = TRUE) {
    x <- as.matrix(x)
    y <- as.numeric(y)
    ynames <- if (is.matrix(y)) {
      rownames(y)
    } else {
      names(y)
    }
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) weights <- rep.int(1, nobs)
    if (is.null(offset)) offset <- rep.int(0, nobs)
    eval(family$initialize)
    # if (coefprior$family == "BIC") coefprior$hyper = as.numeric(nobs)

    newfit <- .Call(C_glm_fit,
      RX = x, RY = y,
      family = family, Roffset = offset,
      Rweights = weights,
      Rpriorcoef = coefprior, Rcontrol = control
    )

    return(newfit)
  }
