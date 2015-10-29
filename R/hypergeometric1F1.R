hypergeometric1F1 = function(a,b,c, laplace=FALSE, log=TRUE) {

    n = length(a);
    out = rep(0, n);
    ans = .C("hypergeometric1F1", as.numeric(a), as.numeric(b), as.numeric(c), out=as.numeric(out), as.integer(n), as.integer(rep(laplace, n)), PACKAGE="BAS")$out
    if (!log) ans = exp(ans)
    return(ans)
}
