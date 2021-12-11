#' Confluent hypergeometric1F1 function
#'
#' Compute the Confluent Hypergeometric function: 1F1(a,b,c,t) =
#' Gamma(b)/(Gamma(b-a)Gamma(a)) Int_0^1 t^(b-1) (1 - t)^(b-a-1) exp(c t) dt
#'
#'
#' @param a arbitrary
#' @param b Must be greater 0
#' @param c arbitrary
#' @param laplace The default is to use the Cephes library; for large a or s
#' this may return an NA, Inf or negative values,, in which case you should use
#' the Laplace approximation.
#' @param log if TRUE, return log(1F1)
#' @author Merlise Clyde (\email{clyde@@stat.duke.edu})
#' @references Cephes library hyp1f1.c
#' @keywords math
#' @examples
#' hypergeometric1F1(11.14756, 0.5, 0.00175097)
#'
#'
#' @rdname hypergeometric1F1
#' @family special functions
#' @export
hypergeometric1F1 = function(a,b,c, laplace=FALSE, log=TRUE) {

    n = length(a);
    out = rep(0, n);
    ans = .C(C_hypergeometric1F1, as.numeric(a), as.numeric(b), as.numeric(c), out=as.numeric(out), as.integer(n),
             as.integer(rep(laplace, n)))$out
    if (!log) ans = exp(ans)
    return(ans)
}
