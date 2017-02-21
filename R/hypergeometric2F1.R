hypergeometric2F1 = function(a,b,c,z, method="Cephes", log=TRUE) {
  out = 1.0
  if (c < b  | b < 0) {
    warning("Must have c > b > 0 in 2F1 function for integral to converge")
    return(Inf)}
  if (abs(z) > 1) {
    warning("integral in 2F1 diverges")
    return(Inf)}
  if (z == 1.0  & c - b - a <= 0 ) ans = Inf
  else {
      if (method == "Laplace") {
        ans = .C(C_logHyperGauss2F1, as.numeric(a), as.numeric(b), as.numeric(c), as.numeric(z), out=as.numeric(out))$out 
        if (!log) ans= exp(ans)
      }
      else {
        ans = .C(C_hypergeometric2F1, as.numeric(a), as.numeric(b), as.numeric(c), as.numeric(z), out=as.numeric(out))$out
        if (is.na(ans)) {
          warning("Cephes routine returned NaN; try Laplace approximation")
          return(NA)}
        else {
          if (ans < 0) 
            warning("2F1 from Cephes library is negative; try Laplace approximation")      
          if (log) ans = log(ans)
        }
      }
    }
      return(ans)
}
