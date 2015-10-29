phi1 = function(a, b, c, x, y) {
# phi1 = int u^{t-1} (1 - v u)^{q - 1} e^-{s u} /B(t, q) (theta
    n = length(t)
    out = rep(0, n);
    ans = .C("phi1", as.numeric(a), as.numeric(b), as.numeric(c), as.numeric(x), as.numeric(y), out=as.numeric(out), as.integer(n), PACKAGE="BAS")$out
    return(ans)
}
