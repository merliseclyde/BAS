#' Gaussian hypergeometric2F1 function
#' 
#' Compute the Gaussian Hypergeometric2F1 function: 2F1(a,b,c,z) = Gamma(b-c)
#' Int_0^1 t^(b-1) (1 - t)^(c -b -1) (1 - t z)^(-a) dt
#' 
#' The default is to use the routine hyp2f1.c from the Cephes library.  If that
#' return a negative value or Inf, one should try method="Laplace" which is
#' based on the Laplace approximation as described in Liang et al JASA 2008.
#' This is used in the hyper-g prior to calculate marginal likelihoods.
#' 
#' @param a arbitrary
#' @param b Must be greater 0
#' @param c Must be greater than b if |z| < 1, and c > b + a if z = 1
#' @param z |z| <= 1
#' @param method The default is to use the Cephes library routine.  This
#' sometimes is unstable for large a or z near one returning Inf or negative
#' values.  In this case, try method="Laplace", which use a Laplace
#' approximation for tau = exp(t/(1-t)).
#' @param log if TRUE, return log(2F1)
#' @return if log=T returns the log of the 2F1 function; otherwise the 2F1
#' function.
#' @author Merlise Clyde (\email{clyde@@stat.duke.edu})
#' @references Cephes library hyp2f1.c
#' 
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O. (2005) Mixtures
#' of g-priors for Bayesian Variable Selection.  Journal of the American
#' Statistical Association.  103:410-423.  \cr
#' \url{http://dx.doi.org/10.1198/016214507000001337}
#' @keywords math
#' @examples
#' hypergeometric2F1(12,1,2,.65)
#' 
#' @rdname hypergeometric2F1
#' @family special functions
#' @export
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
