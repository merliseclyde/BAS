eplogprob = function(lm.obj, thresh=.5, max = 0.99, int=TRUE) {
   pval = summary(lm.obj)$coefficients[,4]
   prob = 1/(1 - exp(1)*pval*log(pval))
   prob[pval > 1/exp(1)] = thresh
   prob[prob > max] = max
   if (int) prob[1] = 1.0
return(prob)
 }

