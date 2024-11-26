# Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
# This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
# License text is available at https://www.gnu.org/licenses/gpl-3.0.html
#

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

    return(parents)
  }


# model.matrix(mf, data)
# attr( , "assign") has where terms are locay ~ ted

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
#' probability, to models that do not include their respective
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
    data = eval(object$call$data, parent.frame())
  )
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

  return(object)
}




# data(Hald)
# bas.hald = bas.lm(Y ~ .^2, data=Hald)
# hald.models = which.matrix(bas.hald$which, n.vars=bas.hald$n.vars)

# par.Hald = make.parents.of.interactions(Y ~ .^2, data=Hald)
# prior = apply(hald.models, 1,
#              FUN=function(x) {prob.hereditary(model=x, parents=par.Hald$parents)})

# .prob.heredity(hald.models[1,], par.Hald$parents)
# force_heredity.bas(bas.hald)

# basd on Don van den Bergh's code for function analytic
# Dedekind numbers, see https://oeis.org/A014466

count.heredity.models <- function(f, max_models = NULL) {75
  t           <- terms(f)
  all_factors <- attr(t, "factors")[-1, , drop = FALSE]
  o           <- attr(t, "order")
  
  # base case 1: all main effects
  if (all(o == 1)) {
    if (!is.null(max_models))
      return(min(2^length(o), max_models))
    else
      return(2^length(o))
  }
  else {
    
    adj_mat <- crossprod(all_factors != 0)
    for (i in seq_len(nrow(adj_mat))) {
      adj_mat[i, ] <- adj_mat[i, ] == o
    }
    tb_o <- tabulate(o)
    
    # m contains the number of models for each order
    m    <- numeric(length(tb_o))
    m[1] <- 2^tb_o[1] # no. models with only main effects
    
    
    Total_Num_Models = m[1]
    # loop upward over higher order interactions
    for (i in 2:length(tb_o)) {
      
      no_lower_order <- sum(tb_o[seq_len(i - 1)])
      
      for (j in seq_len(tb_o[i])) {
        
        # enumerate all combinations of interactions
        combs <- utils::combn(tb_o[i], j)
        for (k in seq_len(ncol(combs))) {
          
          # the subset of the adjacency matrix that corresponds to the current combination
          subset <- adj_mat[no_lower_order + combs[, k], 1:no_lower_order, drop = FALSE]
          
          # all edges imply a parent which must be included (i.e., a constraint)
          constraints <- .colSums(subset, m = j, n = no_lower_order) > 0
          
          # if there are missing constraints, we need to check if they involve
          # interactions of various levels. If they do we need to recurse
          missing <- which(constraints == 0)
          if (any(o[missing] > 1) && length(unique(o[missing])) > 1) {
            
            # recursive case : construct the formula corresponding to the constraint and
            # call analytic again
            # this can probably be sped up by only passing (subsets) of the adjacency matrix around
            cnms <- colnames(adj_mat)[1:no_lower_order]
            new_f_str <- paste("y ~", paste(cnms[missing], collapse = " + "))
            # print(paste(deparse(f), "->", new_f_str))
            new_f <- as.formula(paste("y ~", paste(cnms[missing], collapse = " + ")))
            no_models <- Recall(new_f)
            
          } else {
            
            # base case 2: all free variables are effects of the same order (e.g,. all main effects, all 2-way interactions, etc.)
            no_constraints <- sum(constraints)
            no_models <- 2^(no_lower_order - no_constraints)
            
          }
          
          m[i] <- m[i] + no_models
         
        }
      }
      Total_Num_Models = Total_Num_Models + m[i]
      if (!is.null(max_models) && Total_Num_Models > max_models) {
        return(max_models)
      }
    }
    return(sum(m))
  }
}
