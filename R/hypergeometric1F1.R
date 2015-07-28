hypergeometric1F1 = function(a,b,z, method="Cephes", log=TRUE) {
  out = 1.0
  ans = .C("hypergeometric1F1", as.numeric(a), as.numeric(b), as.numeric(z), out=as.numeric(out, as.integer(length(a)), PACKAGE="BAS")$out
}
  return(ans)
}
