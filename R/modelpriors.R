uniform = function() {
  structure(list(family="Uniform",hyper.parameters=.5), class="prior")
}

Bernoulli=function(probs=0.5) {
  if (length(probs) == 1) {
    if (probs == .5) structure(list(family="Uniform",hyper.parameters=.5), class="prior")
    else structure(list(family="Bernoulli", hyper.parameters=probs), class="prior")}
  else   structure(list(family="Bernoulli", hyper.parameters=probs), class="prior")
}

beta.binomial=function(alpha=1.0, beta=1.0) {
    structure(list(family="Beta-Binomial", hyper.parameters=c(alpha, beta)),
              class="prior")
 }


tr.beta.binomial=function(alpha=1.0, beta=1.0, trunc) {
    structure(list(family="Trunc-Beta-Binomial", hyper.parameters=c(alpha, beta, trunc)),
              class="prior")
 }

tr.power.prior=function(kappa=2, trunc) {
    structure(list(family="Trunc-Power-Prior", hyper.parameters=c(kappa, trunc)),
              class="prior")
 }

tr.poisson=function(lambda, trunc) {
    structure(list(family="Trunc-Poisson", hyper.parameters=c(lambda, trunc)),
              class="prior")
 }
