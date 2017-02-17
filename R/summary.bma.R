print.bas = function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n Marginal Posterior Inclusion Probabilities: \n")
    out = x$probne0
    names(out) = x$namesx
    print.default(format(out, digits = digits), print.gap = 2, 
        quote = FALSE, ...)
    invisible()
}

summary.bas = function(object, n.models = 5, ...) {
  best = order(-object$postprobs)
  n.models = min(n.models, length(best))
  best = best[1:n.models]  
  x = cbind(list2matrix.which(object, best),
               exp(object$logmarg[best] - max(object$logmarg[best])),
               round(object$postprobs[best], 4),
               round(object$R2[best], 4),
               object$size[best],
               object$logmarg[best])
  x = t(x)
  x = cbind(NA, x)
  
  x[1:object$n.vars, 1] = object$probne0
  colnames(x) = c("P(B != 0 | Y)", paste("model", 1:n.models))
  rownames(x) = c(object$namesx, "BF", "PostProbs", "R2", "dim", "logmarg")
  return(x)
  }


