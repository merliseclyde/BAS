#EB.global.EM is now called EB.global

EB.global.bma = function(object, tol= .1, g.0=NULL, max.iterations=100) {
n  = object$n
SSY = var(object$Y)*(n-1)
SSE <- (1.0 - object$R2)*SSY
SSR = SSY - SSE
p = object$size - 1
R2 = object$R2
prior = object$priorprobs

  postmodelprob <- function(R2, p, n, g, prior=1) {
     logmarg  <-  .5*((n - 1 - p)*log(1 + g) - (n-1)*log( 1 + g*(1 - R2)))
     logmarg[p == 0] = 0
     modelprob <- exp(logmarg -  max(logmarg))
     modelprob <- modelprob*prior/sum(modelprob*prior)
     if  (any(is.na(modelprob))) warning("NA's in modelprobs")
    return(modelprob)
   }

best = sort.list(-object$logmarg)[1]
sbest = min(object$shrinkage[best], .99)

if (is.null(g.0)) g.0 = sbest/(1 - sbest)
  tau.0 = g.0 + 1
  phi = (n - 1)/(SSY - (g.0/(1 + g.0))*SSR)
  post.prob = postmodelprob(object$R2,p, n, max(tau.0 - 1, 0), prior)
  tau.0 = sum(post.prob*SSR*phi)/(sum(post.prob*p))
  tau = tau.0 - 2*tol
  it = 0

while (abs(tau - tau.0) > tol | it < max.iterations) {
  g = max(tau - 1, 0)
  phi = (n - 1)/(SSY - (g/(1 + g))*SSR)
  post.prob = postmodelprob(object$R2,p, n, g, prior)
  tau.0 = tau
  tau = sum(post.prob*SSR*phi)/(sum(post.prob*p))
  it = it + 1  
}

g = max(tau -1, 0)
logmarg  =  .5*((n -1 - p)*log(1 + g) - (n-1)*log( 1 + g*(1 - R2)))
logmarg[p == 0] = 0
postprobs = postmodelprob(object$R2,p, n, g, prior)
which = which.matrix(object$which, object$n.var)
object$probne0 = postprobs %*% which 
object$postprobs=postprobs
object$g = g
object$logmarg = logmarg
object$shrinkage = object$shrinkage*0 + g/(1 + g)
object$method = "EB-global"
return(object)
}

