hypergeometric1F1 = function(a,b,z, method="Cephes", log=TRUE) {
  out = 1.0
  if (a < b  | b < 0) {
    warning("Must have a > and b > 0 in 1F1 function for integral to converge")
    return(Inf)}
else{
    ans = .C("hypergeometric1F1", as.numeric(a), as.numeric(b), as.numeric(z), out=as.numeric(out), PACKAGE="BAS")$out
}
  return(ans)
}
