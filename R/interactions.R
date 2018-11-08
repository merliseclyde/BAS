

make.parents.of.interactions <-
  function(mf, data) {
    modelterms <- terms(mf, data = data)
    termnamesX <- attr(modelterms, "term.labels")
    p <- length(termnamesX)
    interactions <- grep(":", termnamesX)
    parents <- diag(p)
    colnames(parents) <- termnamesX
    rownames(parents) <- termnamesX
    for (j in interactions) {
      term <- interactions[j]
      main <- unlist(strsplit(termnamesX[j], ":",
        fixed = TRUE
      ))
      parents.of.term <- main
      for (i in 2:length(main)) {
        parents.of.term <- c(
          parents.of.term,
          utils::combn(main, i, FUN = paste0, collapse = ":")
        )
      }
      parents[j, parents.of.term] <- 1
    }


    X <- model.matrix(modelterms, data = data)

    loc <- attr(X, "assign")[-1] # drop intercept

    parents <- parents[loc, loc]
    parents <- rbind(0, parents)
    parents <- cbind(0, parents)
    parents[1, 1] <- 1
    rownames(parents) <- colnames(X)
    colnames(parents) <- colnames(X)

    # now determine the number of models
    # find the number of terms without any children
    #

    n.childless = sum(colSums(parents) == rowSums(parents))
    term.order = table(lengths(regmatches(termnamesX, gregexpr(":", termnamesX))))
    # colculate the number of models for  childless parents
    n.models = 2^(n.childless - 1)   # drop intercept
    # reduce the number of terms
    n.parents = term.order[1]
    order.of.terms = as.integer(names(term.order))[-1] + 1  # order of interaction

    if (n.parents > 0) {
       model.size = 0:n.parents
       for (i in order.of.terms) {
         n.models = n.models*2^choose(model.size, i)
       }
       n.models = sum(choose(n.parents, model.size)*n.models)
    }
    #
    return(list(parents=parents, n.models = n.models))
  }


# model.matrix(mf, data)
# attr( , "assign") has where terms are located

# mp = BAS:::make.parents.of.interactions(mf, df)



prob.heredity <- function(model, parents, prob = .5) {
  got.parents <- apply(parents, 1,
    FUN = function(x) {
      all(as.logical(model[as.logical(x)]))
    }
  )
  model.prob <- 0
  if (all(model == got.parents)) {
    model.prob <- exp(
      sum(model * log(prob) + (1 - model) * log(1.0 - prob))
    )
  }
  return(model.prob)
}

check.heredity <- function(model, parents, prob = .5) {
  #  p = length(model)  # model has no intercept, while parents does
  #  parents = parents[2:p, 2:p]
  got.parents <- apply(parents, 1,
    FUN = function(x) {
      all(as.logical(model[as.logical(x)]))
    }
  )
  #  browser()
  all(model == got.parents)
}

#' Post processing function to force constraints on interaction inclusion bas BMA objects
#'
#' This function takes the output of a bas object and allows higher
#' order interactions to be included only if their parent
#' lower order interactions terms are in the model, by
#' assigning zero prior probability, and hence posterior
#' probability, to models that do include their respective
#' parents.
#'
#' @param object a bas linear model or generalized linear model object
#' @param prior.prob  prior probability that a term is included conditional on parents being included
#' @return a bas object with updated models, coefficients and summaries obtained removing all models with   zero prior and posterior probabilities.
#' @note Currently prior probabilities are computed using conditional Bernoulli distributions, i.e.  P(gamma_j = 1 | Parents(gamma_j) = 1) = prior.prob.  This is not very efficient for models with a large number of levels.  Future updates will force this at the time of sampling.
#' @author Merlise A Clyde
#' @keywords regression
#' @examples
#'
#' data("chickwts")
#' bas.chk <- bas.lm(weight ~ feed, data = chickwts)
#' #  summary(bas.chk)  # 2^5 = 32 models
#' bas.chk.int <- force.heredity.bas(bas.chk)
#' #  summary(bas.chk.int)  # two models now
#'
#'
#' data(Hald)
#' bas.hald <- bas.lm(Y ~ .^2, data = Hald)
#' bas.hald.int <- force.heredity.bas(bas.hald)
#' image(bas.hald.int)
#'
#' image(bas.hald.int)
#'
#' # two-way interactions
#' data(ToothGrowth)
#' ToothGrowth$dose <- factor(ToothGrowth$dose)
#' levels(ToothGrowth$dose) <- c("Low", "Medium", "High")
#' TG.bas <- bas.lm(len ~ supp * dose, data = ToothGrowth, modelprior = uniform())
#' TG.bas.int <- force.heredity.bas(TG.bas)
#' image(TG.bas.int)
#' @family bas methods
#' @export

force.heredity.bas <- function(object, prior.prob = .5) {
  parents <- make.parents.of.interactions(
             mf = eval(object$call$formula, parent.frame()),
             data = eval(object$call$data, parent.frame()))$parents
  which <- which.matrix(object$which, object$n.vars)
  keep <- apply(which, 1,
    FUN = function(x) {
      check.heredity(model = x, parents = parents)
    }
  )
  #    priorprobs = apply(which, 1,
  #                 FUN=function(x) {prob.heredity(model=x, parents=parents)}
  #    )
  #    keep = priorprobs > 0.0
  object$n.models <- sum(keep)

  object$sampleprobs <- object$sampleprobs[keep] # if method=MCMC ??? reweight
  object$which <- object$which[keep]
  object$priorprobs <- object$priorprobs[keep] / sum(object$priorprobs[keep])
  #    wts = priorprobs[keep]/object$priorprobs[keep]  #importance weights
  wts <- 1
  method <- object$call$method
  if (!is.null(method)) {
    if (method == "MCMC" || method == "MCMC_new") {
      object$freq <- object$freq[keep]
      #      object$postprobs.MCMC = object$freq[keep]*wts
      object$postprobs.MCMC <- object$freq[keep]
      object$postprobs.MCMC <- object$postprobs.MCMC / sum(object$postprobs.MCMC)
      object$probne0.MCMC <- as.vector(object$postprobs.MCMC %*% which[keep, ])
    }
  }
  object$logmarg <- object$logmarg[keep]
  object$shrinkage <- object$shrinkage[keep]
  postprobs.RN <- exp(object$logmarg - min(object$logmarg)) * object$priorprobs
  object$postprobs.RN <- postprobs.RN / sum(postprobs.RN)
  #  browser()
  object$probne0.RN <- as.vector(object$postprobs.RN %*% which[keep, ])

  object$postprobs <- object$postprobs[keep] * wts / sum(object$postprobs[keep] * wts)
  object$probne0 <- as.vector(object$postprobs %*% which[keep, ])

  object$mle <- object$mle[keep]
  object$mle.se <- object$mle.se[keep]
  object$mse <- object$mse[keep]
  object$size <- object$size[keep]
  object$R2 <- object$R2[keep]
  object$df <- object$df[keep]
  object$rank <- object$rank[keep]

  return(object)
}




# data(Hald)
# bas.hald = bas.lm(Y ~ .^2, data=Hald)
# hald.models = which.matrix(bas.hald$which, n.vars=bas.hald$n.vars)

# par.Hald = BAS:::make.parents.of.interactions(Y ~ .^2, data=Hald)
# prior = apply(hald.models, 1,
#              FUN=function(x) {prob.hereditary(model=x, parents=par.Hald$parents)})

# .prob.heredity(hald.models[1,], par.Hald$parents)
# force_heredity.bas(bas.hald)
