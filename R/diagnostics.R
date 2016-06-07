diagnostics = function(obj, type=c("pip","model"),...) {
    if (obj$call$method == "MCMC") {
    if (type == "pip")  {   
    plot(obj$probne0, obj$probs.MCMC,
         xlab="pip (renormalized)",
         ylab="pip (MCMC)", xlim=c(0,1), ylim=c(0,1),
         ...)
    abline(0,1) }
    else {
        plot(obj$postprob, obj$freq/sum(obj$freq),
         xlab="pip (renormalized)",
         ylab="pip (MCMC)",
             ...)
        abline(0,1)
    }
}
    else {
        warning("diagnostic plot is only availble using method='MCMC' for sampling with bas objects.")
    }
}
