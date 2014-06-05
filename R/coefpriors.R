#localEB.g.prior = function() {
# structure(list(family="EB-local-gprior", class="gprior", hyper=NULL), class="prior")
#}

#globalEB.g.prior = function() {
# structure(list(family="EB-global-gprior", class="gprior", hyper=NULL), class="prior")
#}

g.prior = function(g) {
 structure(list(family="fixed-g-prior", g = g, class="gprior", hyper=as.numeric(g)), class="prior")
}


bic.prior = function(n=NULL) {
 structure(list(family="BIC", class="IC", hyper=log(n)), class="prior")
}

aic.prior = function() {
 structure(list(family="AIC", class="IC", hyper=2.0), class="prior")
}

IC.prior = function(penalty) {
 structure(list(family="IC", class="IC", hyper=as.numeric(penalty)), class="prior")
}

