# Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
# This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
# License text is available at https://www.gnu.org/licenses/gpl-3.0.html
#
#' Prediction Method for an Object of Class basglm
#' @description Predictions under model averaging from a BMA (BAS) object for GLMs
#' under different loss functions.
#' @aliases predict.basglm
#' @param object An object of class "basglm", created by \code{bas.glm}
#' @param newdata dataframe, new matrix or vector of data for predictions. May
#' include a column for the intercept or just the predictor variables.  If a
#' dataframe, the variables are extracted using model.matrix using the call
#' that created 'object'.  May be missing in which case the data used for
#' fitting will be used for prediction.
#' @param se.fit indicator for whether to compute se of fitted and predicted
#' values
#' @param type Type of predictions required. The default is  "response" is on the scale of the
#' response variable, with the alternative being on the linear predictor
#' scale, `type ='link'`. Thus for a default binomial model
#' `type = 'response'` gives
#' the predicted probabilities, while with `'link'`, the estimates
#' are of log-odds (probabilities on logit scale).
#' @param top A scalar integer M.  If supplied, calculate results using the subset of the top M models
#' based on posterior probabilities.
#' @param estimator estimator used for predictions.  Currently supported
#' options include: \cr 'HPM' the highest probability model \cr 'BMA' Bayesian
#' model averaging, using optionally only the 'top' models \cr 'MPM' the median
#' probability model of Barbieri and Berger. \cr 'BPM' the model that is
#' closest to BMA predictions under squared error loss. BMA may be computed
#' using only the 'top' models if supplied
#' @param na.action  function determining what should be done with missing values in newdata.
#' The default is to predict NA.
#' @param ... optional extra arguments
#' @return a list of
#'  \item{fit}{predictions using BMA or other estimators}
#'  \item{Ypred}{matrix of predictions under model(s)}
#'  \item{postprobs}{renormalized probabilities of
#' the top models}
#' \item{best}{index of top models included}
#' @details  This function first calls the predict method for class bas
#' (linear models) to form predictions on the linear predictor
#' scale for `BMA`, `HPM`, `MPM` etc. If the estimator is `BMA`
#' and `type='response'` then the
#' inverse link is applied to fitted values for type equal `'link'`
#' and model averaging takes place in the `response` scale. Thus applying
#' the inverse link to BMA estimate with `type = 'link'` is
#' not equal to the fitted values for `type = 'response'` under
#' BMA due to the  nonlinear transformation under the inverse link.
#'
#' @author Merlise Clyde
#' @seealso \code{\link{bas.glm}}, \code{\link{predict.bas}},
#' \code{\link{fitted.bas}}
#' @keywords regression
#' @examples
#'
#'
#' data(Pima.tr, package="MASS")
#' data(Pima.te, package="MASS")
#' Pima.bas = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7, method="BAS",
#'            betaprior=CCH(a=1, b=nrow(Pima.tr)/2, s=0), family=binomial(),
#'            modelprior=uniform())
#' pred = predict(Pima.bas, newdata=Pima.te, top=1)  # Highest Probability model
#' cv.summary.bas(pred$fit, Pima.te$type, score="miss-class")
#'
#' @rdname predict.basglm
#' @family predict methods
#' @family bas methods
#' @export
predict.basglm <- function(object,
                           newdata,
                           se.fit = FALSE,
                           type = c("response", "link"),
                           top = NULL,
                           estimator = "BMA",
                           na.action = na.pass,
                           ...) {
  #    browser()
  if (estimator == "HPM") {
    top <- 1
  }

  # get predictions on linear predictor scale
  pred <- predict.bas(
    object,
    newdata = newdata,
    se.fit = se.fit,
    top = top,
    estimator = estimator,
    na.action = na.action,
    ...
  )

  if (length(type) > 1) {
    type <- type[1]
  }
  #
  # if type is 'link' do not need to do  anything; just return
  # pred at end
  #
  if (type == "response") {
    model.specs <- attributes(pred$fit)
    if (estimator == "BMA") {
      Ypred <- apply(
        pred$Ypred,
        1,
        FUN = function(x) {
          eval(object$family)$linkinv(x)
        }
      )
      if (length(pred$postprobs) > 1) {
        fit <- as.vector(Ypred %*% pred$postprobs)
      } else {
        fit <- as.vector(Ypred)
      }
    }
    else {
      fit <- eval(object$family)$linkinv(pred$fit)
    }
    attributes(fit) <- model.specs

    # standard errors on the response scale using the delta method;
    # mu.eta must be evaluated at the linear predictor (issue #52)
    if (se.fit) {
      # keep results on the scale of the linear predictor so that
      # confint.pred.bas can transform intervals from the link scale
      pred$link <- list(
        fit = as.vector(pred$fit),
        Ypred = pred$Ypred,
        se.fit = pred$se.fit,
        se.pred = pred$se.pred,
        linkinv = eval(object$family)$linkinv
      )
      mu.eta <- eval(object$family)$mu.eta
      if (estimator == "BMA") {
        # rows are models, columns are observations
        dmu <- abs(mu.eta(pred$Ypred))
        pred$se.fit <- pred$se.fit * dmu
        pred$se.pred <- pred$se.pred * dmu
        # variance of the model specific means about the BMA mean
        var.efit <- as.vector(((Ypred - fit)^2) %*% pred$postprobs)
        pred$se.bma.fit <- sqrt(as.vector(pred$postprobs %*% pred$se.fit^2) +
          var.efit)
        pred$se.bma.pred <- sqrt(as.vector(pred$postprobs %*% pred$se.pred^2) +
          var.efit)
      }
      else {
        dmu <- abs(mu.eta(as.vector(pred$fit)))
        pred$se.fit <- pred$se.fit * dmu
        pred$se.pred <- pred$se.pred * dmu
      }
    }

    # replace predictions
    #
    pred$fit <- fit
  }
  return(pred)
}



#' Prediction Method for an object of class BAS
#'
#' Predictions under model averaging or other estimators from a BMA object of
#' class inheriting from 'bas'.
#'
#' Use BMA and/or model selection to form predictions using the top highest
#' probability models.
#'
#' @aliases predict.bas predict
#' @param object An object of class BAS, created by \code{bas}
#' @param newdata dataframe for predictions. If missing, then use the dataframe
#' used for fitting for obtaining fitted and predicted values.
#' @param se.fit indicator for whether to compute se of fitted and predicted
#' values
#' @param type Type of predictions required. The defaults "reponse" is on the scale of
#' the response is the only option currently for linear models (for Gaussian models this
#' is equivalent to type="link").
#' @param top a scalar integer M.  If supplied, subset the top M models, based
#' on posterior probabilities for model predictions and BMA.
#' @param estimator estimator used for predictions.  Currently supported
#' options include: \cr 'HPM' the highest probability model \cr 'BMA' Bayesian
#' model averaging, using optionally only the 'top' models \cr 'MPM' the median
#' probability model of Barbieri and Berger. \cr 'BPM' the model that is
#' closest to BMA predictions under squared error loss. BMA may be computed
#' using only the 'top' models if supplied
#' @param na.action  function determining what should be done with missing values in newdata.
#' The default is to predict NA.
#' @param ... optional extra arguments
#' @return a list of
#' \item{fit}{fitted values based on the selected estimator}
#' \item{Ybma}{predictions using BMA, the same as fit for non-BMA methods for
#' compatibility; will be deprecated}
#' \item{Ypred}{matrix of predictions under
#' each model for BMA}
#' \item{se.fit}{se of fitted values; in the case of BMA
#' this will be a matrix}
#'  \item{se.pred}{se for predicted values; in the case
#' of BMA this will be a matrix}
#'  \item{se.bma.fit}{vector of posterior sd under
#' BMA for posterior mean of the regression function.
#' This will be NULL if estimator is not 'BMA'}
#' \item{se.bma.pred}{vector of posterior sd under BMA
#' for posterior predictive values.  this will be NULL if estimator is not
#' 'BMA'}
#'  \item{best}{index of top models included}
#'  \item{bestmodels}{subset of
#' bestmodels used for fitting or prediction}
#' \item{best.vars}{names of variables in the top model; NULL if estimator='BMA'}
#' \item{df}{scalar or vector of
#' degrees of freedom for models}
#'  \item{estimator}{estimator upon which 'fit'
#' is based.}
#' @author Merlise Clyde
#' @seealso \code{\link{bas}}, \code{\link{fitted.bas}},
#' \code{\link{confint.pred.bas}},  \code{\link{variable.names.pred.bas}}
#' @keywords regression
#' @examples
#'
#' data("Hald")
#' hald.gprior =  bas.lm(Y ~ ., data=Hald, alpha=13, prior="g-prior")
#'
#' predict(hald.gprior, newdata=Hald, estimator="BPM", se.fit=TRUE)
#' # same as fitted
#' fitted(hald.gprior,estimator="BPM")
#' # default is BMA and estimation of mean vector
#' hald.bma = predict(hald.gprior, top=5, se.fit=TRUE)
#' confint(hald.bma)
#'
#' hald.bpm = predict(hald.gprior, newdata=Hald[1,],
#'                     se.fit=TRUE,
#'                     estimator="BPM")
#' confint(hald.bpm)
#' # extract variables
#' variable.names(hald.bpm)
#'
#' hald.hpm = predict(hald.gprior, newdata=Hald[1,],
#'                     se.fit=TRUE,
#'                     estimator="HPM")
#' confint(hald.hpm)
#' variable.names(hald.hpm)
#'
#' hald.mpm = predict(hald.gprior, newdata=Hald[1,],
#'                     se.fit=TRUE,
#'                     estimator="MPM")
#' confint(hald.mpm)
#' variable.names(hald.mpm)
#'
#' @rdname predict.bas
#' @family predict methods
#' @family bas methods
#' @export
predict.bas <- function(object,
                        newdata,
                        se.fit = FALSE,
                        type = c("response","link"),
                        top = NULL,
                        estimator = "BMA",
                        na.action = na.pass,
                        ...) {
  if (!(estimator %in% c("BMA", "HPM", "MPM", "BPM"))) {
    stop("Estimator must be one of 'BMA', 'BPM', 'HPM', or 'MPM'.")
  }

  tt <- terms(object)

  if (missing(newdata) || is.null(newdata)) {
    newdata <- object$X
    insample <- TRUE
  }
  else {
    if (is.data.frame(newdata)) {
      #  newdata = model.matrix(eval(object$call$formula), newdata)
      Terms <- delete.response(tt)
      m <- model.frame(Terms,
        newdata,
        na.action = na.action,
        xlev = object$xlevels
      )
      newdata <- model.matrix(Terms, m,
        contrasts.arg = object$contrasts
      )

      insample <- FALSE
    }
    else {
      stop("use of newdata as a vector is depricated,
           please supply newdata as a dataframe")
      # if (is.vector(newdata)) newdata=matrix(newdata, nrow=1)
    }
  }

  #  browser()
  n <- nrow(newdata)
  if (ncol(newdata) == object$n.vars) {
    newdata <- newdata[, -1, drop = FALSE]
  }  # drop intercept
  if (ncol(newdata) != (object$n.vars - 1)) {
    stop("Dimension of newdata does not match orginal model")
  }
  if (!is.null(object$mean.x)) {
    newdata <- sweep(newdata, 2, object$mean.x)
  }

  df <- object$df


  if (estimator == "MPM") {
    nvar <- object$n.vars - 1
    bestmodel <- (0:nvar)[object$probne0 > .5]
    newX <- cbind(1, newdata)
    best <- 1
    models <- rep(0, nvar + 1)
    models[bestmodel + 1] <- 1
    if (sum(models) > 1) {
      # weights may be a variable in the data frame (as in lm)
      mpm.data <- eval(object$call$data, parent.frame())
      mpm.weights <- eval(object$call$weights, mpm.data, parent.frame())
      if (is.null(mpm.weights)) {
        object <- bas.lm(
          eval(object$call$formula),
          data = mpm.data,
          n.models = 1,
          alpha = object$g,
          initprobs = object$probne0,
          prior = object$prior,
          modelprior = object$modelprior,
          update = NULL,
          bestmodel = models,
          prob.local = .0
        )
      }
      else {
        object <- bas.lm(
          eval(object$call$formula),
          data = mpm.data,
          weights = mpm.weights,
          n.models = 1,
          alpha = object$g,
          initprobs = object$probne0,
          prior = object$prior,
          modelprior = object$modelprior,
          update = NULL,
          bestmodel = models,
          prob.local = .0
        )
      }
      best <- which.max(object$postprobs)
      fit <-
        as.vector(newX[, object$which[[best]] + 1, drop = FALSE] %*% object$mle[[best]]) * object$shrinkage[[best]]
      fit <- fit + (1 - object$shrinkage[[best]]) * (object$mle[[best]])[1]
      # degrees of freedom from the refitted median probability model
      df <- object$df[best]
    }
    else {
      # median probability model is the null model
      best <- which(object$size == 1)[1]
      fit <- rep(as.numeric(object$mle[[best]])[1], nrow(newX))
      df <- object$df[best]
    }
    models <- bestmodel
    attributes(fit) <- list(
      model = models,
      best = best,
      estimator = estimator
    )

    Ybma <- fit
    Ypred <- NULL
    postprobs <- NULL
    best <- NULL
  }
  else {
    if (estimator == "HPM") {
      top <- 1
    }
    postprobs <- object$postprobs
    best <- order(-postprobs)
    if (!is.null(top)) {
      best <- best[1:top]
    }
    models <- object$which[best]
    beta <- object$mle[best]
    gg <- object$shrinkage[best]
    intercept <- object$intercept[best]
    postprobs <- postprobs[best]
    postprobs <- postprobs / sum(postprobs)
    M <- length(postprobs)
    Ypred <- matrix(0, M, n)
    # lm case
    if (is.null(intercept)) {
      for (i in 1:M) {
        beta.m <- beta[[i]]
        model.m <- models[[i]]
        Ypred[i, ] <-
          (newdata[, model.m[-1], drop = FALSE] %*% beta.m[-1]) * gg[i] + beta.m[1]
      }
    }
    else {
      for (i in 1:M) {
        beta.m <- beta[[i]]
        model.m <- models[[i]]
        Ypred[i, ] <-
          (newdata[, model.m[-1], drop = FALSE] %*% beta.m[-1]) * gg[i] + intercept[i]
      }
    }


    df <- df[best]
    Ybma <- t(Ypred) %*% postprobs
    fit <- as.vector(Ybma)
    if (estimator == "HPM") {
      models <- unlist(object$which[best])
      attributes(fit) <- list(
        model = models,
        best = best,
        estimator = estimator
      )
    }
    if (estimator == "BPM") {
      #    browser()
      dis <- apply(
        sweep(Ypred, 2, Ybma),
        1,
        FUN = function(x) {
          x[is.na(x)] <- 0 # ignore NA's in finding closest model
          sum(x^2)
        }
      )
      bestBPM <- which.min(dis)
      fit <- as.vector(Ypred[bestBPM, ])
      models <- unlist(object$which[best[bestBPM]])
      best <- best[bestBPM]
      # df has already been subset and reordered by best
      df <- df[bestBPM]
      attributes(fit) <- list(
        model = models,
        best = best,
        estimator = estimator
      )
    }
  }
  # browser()
  se <- list(
    se.fit = NULL,
    se.pred = NULL,
    se.bma.fit = NULL,
    se.bma.pred = NULL
  )

  if (se.fit) {
    if (estimator != "BMA") {
      se <- .se.fit(fit, newdata, object, insample)
    }
    else {
      se <- .se.bma(
        Ybma, newdata, Ypred, best, object,
        insample
      )
    }
  }

  best.vars <- object$namesx # BMA case
  if (!is.list(models)) {
    best.vars <- object$namesx[models + 1]
  }


  out <- list(
    fit = fit,
    Ybma = Ybma,
    Ypred = Ypred,
    postprobs = postprobs,
    se.fit = se$se.fit,
    se.pred = se$se.pred,
    se.bma.fit = se$se.bma.fit,
    se.bma.pred = se$se.bma.pred,
    df = df,
    best = best,
    bestmodel = models,
    best.vars = best.vars,
    estimator = estimator
  )

  class(out) <- "pred.bas"
  return(out)
}




#' Fitted values for a BAS BMA objects
#'
#' Calculate fitted values for a BAS BMA object
#'
#' Calculates fitted values at observed design matrix using either the highest
#' probability model, 'HPM', the posterior mean (under BMA) 'BMA', the median
#' probability model 'MPM' or the best predictive model 'BPM".  The median
#' probability model is defined by including variable where the marginal
#' inclusion probability is greater than or equal to 1/2. For type="BMA", the
#' weighted average may be based on using a subset of the highest probability
#' models if an optional argument is given for top.  By default BMA uses all
#' sampled models, which may take a while to compute if the number of variables
#' or number of models is large.  The "BPM" is found be computing the squared
#' distance of the vector of fitted values for a model and the fitted values
#' under BMA and returns the model with the smallest distance.  In the presence
#' of multicollinearity this may be quite different from the MPM, with extreme
#' collinearity may drop relevant predictors.
#'
#' @aliases fitted.bas fitted
#' @param object An object of class 'bas' as created by \code{\link{bas}}
#' @param type type equals "response" or "link" in the case of GLMs (default is 'link')
#' @param estimator estimator type of fitted value to return. Default is to use
#' BMA with all models. Options include \cr 'HPM' the highest probability model
#' \cr 'BMA' Bayesian model averaging, using optionally only the 'top' models
#' \cr 'MPM' the median probability model of Barbieri and Berger.  'BPM' the
#' model that is closest to BMA predictions under squared error loss
#' @param top optional argument specifying that the 'top' models will be used
#' in constructing the BMA prediction, if NULL all models will be used.  If
#' top=1, then this is equivalent to 'HPM'
#' @param na.action function determining what should be done with missing values in newdata. The default is to predict NA.
#' @param ... optional arguments, not used currently
#' @return A vector of length n of fitted values.
#' @author Merlise Clyde \email{clyde@duke.edu}
#' @seealso \code{\link{predict.bas}}  \code{\link{predict.basglm}}
#' @references Barbieri, M.  and Berger, J.O. (2004) Optimal predictive model
#' selection. Annals of Statistics. 32, 870-897. \cr
#' \url{https://projecteuclid.org/euclid.aos/1085408489&url=/UI/1.0/Summarize/euclid.aos/1085408489}
#'
#' Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling for
#' Variable Selection and Model Averaging. Journal of Computational Graphics
#' and Statistics.  20:80-101 \cr
#' \doi{10.1198/jcgs.2010.09049}
#' @keywords regression
#' @examples
#'
#' data(Hald)
#' hald.gprior =  bas.lm(Y~ ., data=Hald, prior="ZS-null", initprobs="Uniform")
#' plot(Hald$Y, fitted(hald.gprior, estimator="HPM"))
#' plot(Hald$Y, fitted(hald.gprior, estimator="BMA", top=3))
#' plot(Hald$Y, fitted(hald.gprior, estimator="MPM"))
#' plot(Hald$Y, fitted(hald.gprior, estimator="BPM"))
#'
#' @rdname fitted
#' @family bas methods
#' @family predict methods
#' @export
fitted.bas <- function(object,
                       type = "link",
                       estimator = "BMA",
                       top = NULL,
                       na.action = na.pass,
                       ...) {
  nmodels <- length(object$which)
  X <- object$X
  if (is.null(top)) {
    top <- nmodels
  }
  if (estimator == "HPM") {
    yhat <- predict(
      object,
      newdata = NULL,
      top = 1,
      estimator = "HPM", type = type,
      na.action = na.action
    )$fit
  }
  if (estimator == "BMA") {
    yhat <- predict(
      object,
      newdata = NULL,
      top = top,
      estimator = "BMA", type = type,
      na.action = na.action
    )$fit
  }
  if (estimator == "MPM") {
    yhat <- predict(
      object,
      newdata = NULL,
      top = top,
      estimator = "MPM", type = type,
      na.action = na.action
    )$fit
  }
  if (estimator == "BPM") {
    yhat <- predict(
      object,
      newdata = NULL,
      top = top,
      estimator = "BPM", type = type,
      na.action = na.action
    )$fit
  }

  return(as.vector(yhat))
}

# Internal helpers for standard errors of fitted values (GitHub issue #52)
#
# The covariance matrix of the coefficients in model M is proportional to
# (X_M^T W X_M)^{-1}, where W holds the prior weights for linear models and
# the prior weights times the IRLS working weights evaluated at the MLE for
# GLMs.  Earlier versions used the unweighted (X_M^T X_M)^{-1} so that
# se.fit and se.pred were incorrect whenever W was not the identity.

# family object for a BAS object (NULL for linear models)
.bas.family <- function(object) {
  if (!inherits(object, "basglm") || is.null(object$family)) {
    return(NULL)
  }
  family <- eval(object$family)
  if (is.function(family)) {
    family <- family()
  }
  family
}

# weights W that define the information matrix X^T W X for a model, and
# the dispersion for GLMs
.info.weights <- function(object, model, best) {
  nobs <- nrow(object$X)
  prior.weights <- object$weights
  if (is.null(prior.weights)) {
    prior.weights <- rep(1, nobs)
  }
  family <- .bas.family(object)
  if (is.null(family)) {
    return(list(weights = prior.weights, dispersion = NULL))
  }

  offset <- object$offset
  if (is.null(offset)) {
    offset <- rep(0, nobs)
  }
  eta <- as.vector(object$X[, model + 1, drop = FALSE] %*% object$mle[[best]]) +
    offset
  mu <- family$linkinv(eta)
  varmu <- family$variance(mu)
  wts <- prior.weights * family$mu.eta(eta)^2 / varmu

  dispersion <- 1
  if (!(family$family %in% c("binomial", "poisson"))) {
    # Pearson estimate as in summary.glm
    pos <- prior.weights > 0
    df.r <- sum(pos) - length(model)
    dispersion <- if (df.r > 0) {
      sum((prior.weights * (object$Y - mu)^2 / varmu)[pos]) / df.r
    } else {
      NaN
    }
  }
  list(weights = wts, dispersion = dispersion)
}

# diagonal of Xnew (X^T W X)^{-1} Xnew^T using the pivoted Cholesky
# decomposition of X^T W X (LAPACK dpstrf, as in cholregpivot in C for
# bas.lm); aliased columns are dropped if X^T W X is not full rank
# (GitHub issue #70)
.quadform.XtWXinv <- function(Xnew, XtWX, tol = 1e-7) {
  # chol warns if the matrix is rank deficient
  R <- suppressWarnings(chol(XtWX, pivot = TRUE, tol = tol))
  r <- attr(R, "rank")
  piv <- attr(R, "pivot")[seq_len(r)]
  R <- R[seq_len(r), seq_len(r), drop = FALSE]
  RinvTx <- backsolve(R, t(Xnew[, piv, drop = FALSE]), transpose = TRUE)
  colSums(RinvTx^2)
}

# scale factor s so that Var(fitted value) = sigma^2 * s for one model
#
# With the predictors centered by their W-weighted means the intercept is
# orthogonal to the other coefficients, with variance sigma^2/sum(W), and
# the posterior covariance of the remaining coefficients is
# shrinkage * sigma^2 (Xc^T W Xc)^{-1}.  Since
#  x^T (X^T W X)^{-1} x = 1/sum(W) + xc^T (Xc^T W Xc)^{-1} xc
# the scale is  1/sum(W) + shrinkage * (x^T (X^T W X)^{-1} x - 1/sum(W))
#
# Xnew and Xold both include the column of ones; the quadratic form does
# not depend on how the remaining columns are centered as long as the same
# centering is applied to both.
#
# For linear models W does not depend on the model, so that X^T W X for the
# full model may be computed once (see .XtWX.lm) and passed as XtWX; the
# cross-products for a model are then a submatrix.  For GLMs W depends on
# the MLE of the model and the cross-products are computed for each model.
.scale.fit <- function(model, best, Xnew, Xold, object, XtWX = NULL) {
  info <- .info.weights(object, model, best)
  sumw <- sum(info$weights)
  cols <- model + 1
  if (is.null(XtWX) || !is.null(info$dispersion)) {
    XtWX.m <- crossprod(sqrt(info$weights) * Xold[, cols, drop = FALSE])
  }
  else {
    XtWX.m <- XtWX[cols, cols, drop = FALSE]
  }
  xiXTXxiT <- .quadform.XtWXinv(Xnew[, cols, drop = FALSE], XtWX.m) -
    1 / sumw
  list(
    scale = 1 / sumw + object$shrinkage[best] * xiXTXxiT,
    dispersion = info$dispersion
  )
}

# X^T W X for the full model in linear models; NULL for GLMs
.XtWX.lm <- function(Xold, object) {
  if (!is.null(.bas.family(object))) {
    return(NULL)
  }
  wts <- object$weights
  if (is.null(wts)) {
    wts <- rep(1, nrow(Xold))
  }
  crossprod(sqrt(wts) * Xold)
}

# design matrix used in fitting, centered in the same way as newdata in
# predict.bas
.Xold.centered <- function(object) {
  Xold <- object$X[, -1, drop = FALSE]
  if (!is.null(object$mean.x)) {
    Xold <- sweep(Xold, 2, object$mean.x)
  }
  cbind(1, Xold)
}

# weighted total sum of squares for linear models
.weighted.ssy <- function(object) {
  wts <- object$weights
  if (is.null(wts)) {
    wts <- rep(1, length(object$Y))
  }
  sum(wts * (object$Y - weighted.mean(object$Y, wts))^2)
}

.se.fit <- function(yhat, X, object, insample) {
  model <- attr(yhat, "model")
  best <- attr(yhat, "best")

  df <- object$df[best]
  shrinkage <- object$shrinkage[best]

  sf <- .scale.fit(model, best, cbind(1, X), .Xold.centered(object), object)
  scale_fit <- sf$scale

  if (is.null(sf$dispersion)) {
    # linear model
    bayes_mse <- .weighted.ssy(object) *
      (1 - shrinkage * object$R2[best]) / df
  }
  else {
    bayes_mse <- sf$dispersion
  }
  # se.pred is for a new observation with unit weight
  se.fit <- sqrt(bayes_mse * scale_fit)
  se.pred <- sqrt(bayes_mse * (1 + scale_fit))
  return(list(
    se.fit = se.fit,
    se.pred = se.pred,
    residual.scale = sqrt(bayes_mse)
  ))
}

.se.bma <- function(fit, Xnew, Ypred, best, object, insample) {
  df <- object$df[best]
  shrinkage <- object$shrinkage[best]

  Xnew <- cbind(1, Xnew)
  Xold <- .Xold.centered(object)
  XtWX <- .XtWX.lm(Xold, object)

  sf <- lapply(seq_along(best), FUN = function(i) {
    .scale.fit(object$which[[best[i]]], best[i], Xnew, Xold, object, XtWX)
  })
  # rows are observations, columns are models
  scale_fit <- matrix(
    unlist(lapply(sf, function(x) x$scale)),
    nrow = nrow(Xnew), ncol = length(best)
  )

  if (is.null(.bas.family(object))) {
    bayes_mse <- .weighted.ssy(object) *
      (1 - shrinkage * object$R2[best]) / df
  }
  else {
    bayes_mse <- unlist(lapply(sf, function(x) x$dispersion))
  }

  var.fit <- sweep(scale_fit, 2, bayes_mse, FUN = "*")
  var.pred <- sweep((1 + scale_fit), 2, bayes_mse, FUN = "*")


  postprobs <- object$postprobs[best]

  # expected variance
  evar.fit <- as.vector(var.fit %*% postprobs)
  evar.pred <- as.vector(var.pred %*% postprobs)
  # variance of expectations
  var.efit <- as.vector(postprobs %*% (sweep(Ypred, 2, fit))^2)

  se.fit <- sqrt(evar.fit + var.efit)
  se.pred <- sqrt(evar.pred + var.efit)


  return(
    list(
      se.bma.fit = se.fit,
      se.bma.pred = se.pred,
      se.fit = t(sqrt(var.fit)),
      se.pred = t(sqrt(var.pred)),
      residual.scale = sqrt(bayes_mse)
    )
  )
}

#' Extract the variable names for a model from a BAS prediction object
#'
#' @description S3 method for class 'pred.bas'.  Simple utility
#' function to extract the variable names.  Used to print names
#' for the selected models using estimators for 'HPM', 'MPM' or 'BPM".
#' for the selected model created by \code{predict} for BAS
#' objects.
#' @param object a BAS object created by \code{predict} from a BAS
#' `bas.lm` or `bas.glm` object
#' @param ...  other arguments to pass on
#' @return a character vector with the names of the variables
#' included in the selected model; in the case of 'BMA' this will
#' be all variables
#' @seealso \code{\link{predict.bas}}
#' @method variable.names pred.bas
#' @rdname variable.names.pred.bas
#' @aliases variable.names.pred.bas variable.names
#' @family predict methods
#' @family bas methods
#' @examples
#' data(Hald)
#' hald.gprior =  bas.lm(Y~ ., data=Hald, prior="ZS-null", modelprior=uniform())
#' hald.bpm = predict(hald.gprior, newdata=Hald[1,],
#'                    se.fit=TRUE,
#'                    estimator="BPM")
#' variable.names(hald.bpm)
#' @export
#'
variable.names.pred.bas <- function(object, ...) {
  if (inherits(object, "pred.bas")) {
    object$best.vars
  }
}
