#' eplogprob - Compute approximate marginal inclusion probabilities from
#' pvalues
#' 
#' \code{eplogprob} calculates approximate marginal posterior inclusion
#' probabilities from p-values computed from a linear model using a lower bound
#' approximation to Bayes factors.  Used to obtain initial inclusion
#' probabilities for sampling using Bayesian Adaptive Sampling \code{bas.lm}
#' 
#' Sellke, Bayarri and Berger (2001) provide a simple calibration of p-values
#' 
#' BF(p) = -e p log(p)
#' 
#' which provide a lower bound to a Bayes factor for comparing H0: beta = 0
#' versus H1: beta not equal to 0, when the p-value p is less than 1/e.  Using
#' equal prior odds on the hypotheses H0 and H1, the approximate marginal
#' posterior inclusion probability
#' 
#' p(beta != 0 | data ) = 1/(1 + BF(p))
#' 
#' When p > 1/e, we set the marginal inclusion probability to 0.5 or the value
#' given by \code{thresh}.
#' 
#' @param lm.obj a linear model object
#' @param thresh the value of the inclusion probability when if the p-value >
#' 1/exp(1), where the lower bound approximation is not valid.
#' @param max maximum value of the inclusion probability; used for the
#' \code{bas.lm} function to keep initial inclusion probabilities away from 1.
#' @param int If the Intercept is included in the linear model, set the
#' marginal inclusion probability corresponding to the intercept to 1
#' @return \code{eplogprob} returns a vector of marginal posterior inclusion
#' probabilities for each of the variables in the linear model. If int = TRUE,
#' then the inclusion probability for the intercept is set to 1.  If the model
#' is not full rank, variables that are linearly dependent base on the QR
#' factorization will have NA for their p-values.  In bas.lm, where the
#' probabilities are used for sampling, the inclusion probability is set to 0.
#' @author Merlise Clyde \email{clyde@@stat.duke.edu}
#' @seealso \code{\link{bas}}
#' @references Sellke, Thomas, Bayarri, M. J., and Berger, James O.  (2001),
#' ``Calibration of p-values for testing precise null hypotheses'', The
#' American Statistician, 55, 62-71.
#' @keywords regression
#' @examples
#' 
#' library(MASS)
#' data(UScrime)
#' UScrime[,-2] = log(UScrime[,-2])
#' eplogprob(lm(y ~ ., data=UScrime))
#' 
eplogprob = function(lm.obj, thresh=.5, max = 0.99, int=TRUE) {
    pval = summary(lm.obj)$coefficients[,4]
    prob = 1/(1 - exp(1)*pval*log(pval))
    prob[pval > 1/exp(1)] = thresh
    prob[prob > max] = max
    if (int) prob[1] = 1.0
    if (any(is.na(prob))) {
        warning("Model is not full rank, do not use eplogp approximation to start sampling\n")
    }                    
return(prob)
}



#' eplogprob.marg - Compute approximate marginal inclusion probabilities from
#' pvalues
#' 
#' \code{eplogprob.marg} calculates approximate marginal posterior inclusion
#' probabilities from p-values computed from a series of simple linear
#' regression models using a lower bound approximation to Bayes factors.  Used
#' to order variables and if appropriate obtain initial inclusion probabilities
#' for sampling using Bayesian Adaptive Sampling \code{bas.lm}
#' 
#' Sellke, Bayarri and Berger (2001) provide a simple calibration of p-values
#' 
#' BF(p) = -e p log(p)
#' 
#' which provide a lower bound to a Bayes factor for comparing H0: beta = 0
#' versus H1: beta not equal to 0, when the p-value p is less than 1/e.  Using
#' equal prior odds on the hypotheses H0 and H1, the approximate marginal
#' posterior inclusion probability
#' 
#' p(beta != 0 | data ) = 1/(1 + BF(p))
#' 
#' When p > 1/e, we set the marginal inclusion probability to 0.5 or the value
#' given by \code{thresh}. For the eplogprob.marg the marginal p-values are
#' obtained using statistics from the p simple linear regressions
#' 
#' P(F > (n-2) R2/(1 - R2)) where F ~ F(1, n-2) where R2 is the square of the
#' correlation coefficient between y and X_j.
#' 
#' @param Y response variable
#' @param X design matrix with a column of ones for the intercept
#' @param thresh the value of the inclusion probability when if the p-value >
#' 1/exp(1), where the lower bound approximation is not valid.
#' @param max maximum value of the inclusion probability; used for the
#' \code{bas.lm} function to keep initial inclusion probabilities away from 1.
#' @param int If the Intercept is included in the linear model, set the
#' marginal inclusion probability corresponding to the intercept to 1
#' @return \code{eplogprob.prob} returns a vector of marginal posterior
#' inclusion probabilities for each of the variables in the linear model. If
#' int = TRUE, then the inclusion probability for the intercept is set to 1.
#' @author Merlise Clyde \email{clyde@@stat.duke.edu}
#' @seealso \code{\link{bas}}
#' @references Sellke, Thomas, Bayarri, M. J., and Berger, James O.  (2001),
#' ``Calibration of p-values for testing precise null hypotheses'', The
#' American Statistician, 55, 62-71.
#' @keywords regression
#' @examples
#' 
#' library(MASS)
#' data(UScrime)
#' UScrime[,-2] = log(UScrime[,-2])
#' eplogprob(lm(y ~ ., data=UScrime))
#' 
eplogprob.marg = function(Y,X, thresh=.5, max = 0.99, int=TRUE) {
#    browser()
    R2 = apply(X[,-1], 2, FUN=function(x, y=Y) {cor(x,y)^2})
    n = length(Y)
    Fstat = (n-2)*R2/(1 - R2)
    pval = 1 - pf(Fstat, 1, n-2)
    prob = 1/(1 - exp(1)*pval*log(pval))
    prob[pval > 1/exp(1)] = thresh
    prob[pval < 10^-16] = max
    prob[prob > max] = max
    if (int) prob = c(1.0, prob)
                 
return(prob)
}


