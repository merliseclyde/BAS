#' Extract the Median Probability Model
#' @description Extracts the Median Probability Model from a bas object
#' @param object An object of class "bas" or "basglm"
#' @return  a new object with of class "bas" or "basglm" with the Median
#' Probability Model
#' @details The Median Probability Model is the model where variables are
#' included if the marginal posterior probabilty of the coefficient being
#' zero is greater than 0.5.  As this model may not have been sampled (and even
#' if it has) it is oftern faster to refit the model using bas, rather than
#' search the list of models to see where it was included.
#' @examples
#' data(Hald, package=BAS)
#' hald_bic =  bas.lm(Y ~ ., data=Hald, alpha=13, prior="BIC")
#' extract_MPM(hald_bic)
#'
#' data(Pima.tr, package="MASS")
#' Pima_bas = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7, method="BAS",
#'            betaprior=CCH(a=1, b=nrow(Pima.tr)/2, s=0), family=binomial(),
#'            modelprior=uniform())
#' extract_MPM(Pima_bas)
#' @family bas methods
#' @export
extract_MPM = function(object) {
#  if (!(class(object) %in% c("basglm", "bas"))) {
#    stop("requires an object of class 'bas' or 'basglm'") }
  nvar <- object$n.vars - 1
  bestmodel <- as.numeric(object$probne0 > .5)

  if (is.null(object$call$weights)) {
      object$call$weights = NULL }

  if ( !("basglm" %in% class(object))) {
    # call lm
      newobject <- bas.lm(
        eval(object$call$formula),
        data = eval(object$call$data, parent.frame()),
        weights = eval(object$call$weights),
        n.models = 1,
        alpha = object$g,
        initprobs = object$probne0,
        prior = object$prior,
        modelprior = object$modelprior,
        update = NULL,
        bestmodel = bestmodel
      )

  }
else {
    glm_family = eval(object$family, parent.frame())$family
    family <- get(glm_family, mode = "function", envir = parent.frame())
    newobject <- bas.glm(
      eval(object$call$formula),
      data = eval(object$call$data, parent.frame()),
      weights =  eval(object$call$weights),
      family = family,
      n.models = 1L,
      initprobs = object$probne0,
      betaprior = object$betaprior,
      modelprior = object$modelprior,
      update = NULL,
      bestmodel = bestmodel
  )
}
  newobject$probne0 = object$probne0
  mf = object$call
  mf$n.models = 1
  mf$bestmodel = bestmodel
  newobject$call = mf
  return(newobject)
}
