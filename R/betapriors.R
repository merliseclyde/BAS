EB.local = function() {
    structure(list(family="EB-local", class="EB", hyper.parameters=list(local=TRUE)),
    class="prior")
}

CCH = function(alpha, beta, s=0) {
#    if (beta == 2 & alpha == 2 & s == 0)   {
#        structure(list(family="Truncated-Gamma", class="TCCH", hyper.parameters=NULL),
#                       class="prior")}
#    else {      
        structure(list(family="CCH", class="TCCH",
                       hyper.parameters=list(alpha=alpha, beta=beta, s=s)),
                  class="prior")
    #}
}

tCCH = function(alpha=1, beta=2, s=0, r=3/2, v=1, theta=1) {
#    if (beta == 2 & alpha == 2 & s == 0)   {
#        structure(list(family="Truncated-Gamma", class="TCCH", hyper.parameters=NULL),
#                       class="prior")}
#    else {      
        structure(list(family="tCCH", class="TCCH",
                       hyper.parameters=list(alpha=alpha, beta=beta, s=s,
                           r=r, v=v, theta=theta)),
                  class="prior")
    #}
    }

intrinsic = function(n) {
#    if (beta == 2 & alpha == 2 & s == 0)   {
#        structure(list(family="Truncated-Gamma", class="TCCH", hyper.parameters=NULL),
#                       class="prior")}
#    else {      
        structure(list(family="intrinsic", class="TCCH",
                       hyper.parameters=list(alpha=1, beta=1, s=0, r=1, n=n)),
                  class="prior")
    #}
    }

hyper.g.n = function(alpha=3, n) {
#    if (beta == 2 & alpha == 2 & s == 0)   {
#        structure(list(family="Truncated-Gamma", class="TCCH", hyper.parameters=NULL),
#                       class="prior")}
#    else {      
        structure(list(family="tCCH", class="TCCH",
                       hyper.parameters=list(alpha=alpha-2, beta=2, s=0,
                           r=alpha/2, v=1, theta=1/n)),
                  class="prior")
    #}
}

Jeffreys = function() {
        structure(list(family="Jeffreys", class="TCCH",
                       hyper.parameters=list(alpha=0, beta=2, s=0.0)),
                  class="prior")
    }


hyper.g = function(alpha=3) {
    if (alpha <= 2 )   {
        return("alpha must be greater than 2 in hyper.g prior")
    }
    
    else {      
        structure(list(family="CCH", class="TCCH",
                       hyper.parameters=list(alpha=alpha-2.0, beta=2, s=0.0)),
                  class="prior")}
}

TG = function(alpha=2) {
    structure(list(family="TG", class="TCCH",
                       hyper.parameters=list(alpha=alpha, beta=2, s=0.0)),
                  class="prior")
}


beta.prime=function(n=NULL) {
    structure(list(family="betaprime", class="TCCH", hyper.parameters=list(n=n, alpha=.5)),
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
    structure(list(family="g.prior", g = as.numeric(g), class="g-prior",
                   hyper=as.numeric(g),
                   hyper.parameters= list(g=g)),
              class="prior")
}
