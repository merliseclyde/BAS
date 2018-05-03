#' BAS MCMC diagnostic plot
#'
#' Function to help assess convergence of MCMC sampling for bas ojects.
#'
#' BAS calculates posterior model probabilites in two ways when method="MCMC".
#' The first is using the relative Monte Carlo frequencies of sampled models.
#' The second is to renormalize the marginal likelihood times prior
#' probabilities over the sampled models.  If the markov chain has converged,
#' these two quantities should be the same and fall on a 1-1 line.  If not,
#' running longer may be required.  If the chain has not converged, the Monte
#' Carlo frequencise may have have less bias, although may exhibit more
#' variability.
#'
#' @param obj an object created by bas.lm or bas.glm
#' @param type type of diagnostic plot.  If "pip" the marginal inclusion
#' probabilities are used, while if "model", plot posterior model probabilities
#' @param ... additional graphics parameters to be passed to plot
#' @return a plot with of the marginal inclusion probabilities (pip) estimated
#' by MCMC and renormalized marginal likelihoods times prior probabilities or
#' model probabilities.
#' @author Merlise Clyde (\email{clyde@duke.edu})
#' @examples
#'
#' library(MASS)
#' data(UScrime)
#' UScrime[,-2] = log(UScrime[,-2])
#' crime.ZS =  bas.lm(y ~ .,
#'                    data=UScrime,
#'                    prior="ZS-null",
#'                    modelprior=uniform(),
#'                    method = "MCMC",
#'                    MCMC.iter = 1000)   # short run for the example
#' diagnostics(crime.ZS)
#'
#' @family bas methods
#' @export
diagnostics = function(obj, type=c("pip","model"),...) {
    if (obj$call$method == "MCMC") {
    for (i in 1:length(type)) {
    if (type[i] == "pip")  {
      plot(obj$probne0.RN, obj$probne0.MCMC,
          xlab="pip (renormalized)",
          ylab="pip (MCMC)", xlim=c(0,1), ylim=c(0,1),
          main="Convergence Plot: Posterior Inclusion Probabilities",
          ...)
      abline(0,1) }
    else {
        ax.lim = range(pretty(c(obj$postprobs.RN, obj$postprobs.MCMC)))
        plot(obj$postprobs.RN, obj$postprobs.MCMC,
         xlab="p(M | Y) (renormalized)",
         ylab="p(M | Y) (MCMC)", xlim=ax.lim, ylim=ax.lim,
         main="Convergence Plot: Posterior Model Probabilities",
             ...)
        abline(0,1)
    }
    }
}
    else {
        simpleError("diagnostic plot is only availble using method='MCMC' for sampling with bas objects.")
    }
}
