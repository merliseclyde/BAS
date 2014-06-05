print.bma = function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n Marginal Posterior Inclusion Probabilities: \n")
    out = x$probne0
    names(out) = x$namesx
    print.default(format(out, digits = digits), print.gap = 2, 
        quote = FALSE, ...)
    invisible()
}

summary.bma = function(object, n.models = 5, ...) {
  best = order(-object$postprobs)
  best = best[1:n.models]  
  x = cbind(list2matrix.which(object, best),
               exp(object$logmarg[best] - max(object$logmarg[best])),
               round(object$postprobs[best], 4),
               round(object$R2[best], 4),
               object$size[best],
               object$logmarg[best])

  colnames(x) = c(object$namesx, "BF", "PostProbs", "R2", "dim", "logmarg")
  x}


