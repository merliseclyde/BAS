hypergeometric1F1 = function(a,b,z, laplace=TRUE, log=TRUE) {

    n = length(a);
    out = rep(0, n);
    ans = .C("hypergeometric1F1", as.numeric(a), as.numeric(b), as.numeric(z), out=as.numeric(out), as.integer(n), as.integer(rep(laplace, n)), PACKAGE="BAS")$out
    return(ans)
}
