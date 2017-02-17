diagnostics = function(obj, type=c("pip","model"),...) {
    if (obj$call$method == "MCMC") {
    for (i in 1:length(type)) {  
    if (type[i] == "pip")  {   
      plot(obj$probne0.RN, obj$probne0.MCMC,
          xlab="pip (renormalized)",
          ylab="pip (MCMC)", xlim=c(0,1), ylim=c(0,1),
          ...)
      abline(0,1) }
    else {
        plot(obj$postprobs.RN, obj$postprobs.MCMC,
         xlab="pip (renormalized)",
         ylab="pip (MCMC)",
             ...)
        abline(0,1)
    }
    }
}
    else {
        warning("diagnostic plot is only availble using method='MCMC' for sampling with bas objects.")
    }
}
