CCH = function(alpha, beta, s=0) {
    if (beta == 2 & alpha == 2 & s == 0)   {
        structure(list(family="Truncated-Gamma", class="TCCH", hyper.parameters=NULL),
                       class="prior")}
    else {      
        structure(list(family="CCH", class="TCCH",
                       hyper.parameters=list(alpha=alpha, beta=beta, s=s)),
                  class="prior")}
}


beta.prime=function(n=NULL) {
    structure(list(family="betaprime", class="TCCH", hyper.parameters=list(n=n, alpha=.25)),
              class="prior")
}

robust = function(n=NULL) {
    structure(list(family="robust", class="TCCH",  hyper.parameters = list(n=n)),
                   class="prior")
}

bic.prior = function(n=NULL) {
    structure(list(family="BIC", class="IC", hyper.parameters = list(penalty=log(n)),
                   hyper=log(n)),
                   class="prior")
}

aic.prior = function() {
    structure(list(family="AIC", class="IC", hyper.parameters = list(penalty=2),
                   hyper=2),
                   class="prior")
}

    
IC.prior = function(penalty) {
    structure(list(family="IC", class="IC", hyper=as.numeric(penalty),
                   hyper.parameters = list(penalty=penalty)),
              class="prior")
}    

g.prior = function(g) {
    structure(list(family="fixed-g-prior", g = g, class="gprior", hyper=as.numeric(g),
                   hyper.parameters= list(g=g)),
              class="prior")
}
