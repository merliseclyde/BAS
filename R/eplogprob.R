#' @export
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


