bayesglm.fit <-
function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, nobs), family = binomial(),
            coefprior = bic.prior(nobs),
            control = glm.control(),intercept=TRUE) 
{

  x <- as.matrix(x)
  y <- as.numeric(y)
  ynames <- if (is.matrix(y))     rownames(y)
            else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights))   weights <- rep.int(1, nobs)
  if (is.null(offset))    offset <- rep.int(0, nobs)
  eval(family$initialize)
  if (coefprior$family == "BIC") coefprior$hyper = as.numeric(nobs)

  newfit = .Call("glm_fit",
    RX=x, RY = y,
    family=family, Roffset = offset,
    Rweights = weights,
      Rpriorcoef = coefprior, Rcontrol=control, PACKAGE="BAS")

  return(newfit)
}
