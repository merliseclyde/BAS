#' Compound Confluent hypergeometric function of two variables
#'
#' Compute the Confluent Hypergeometric function of two variables, also know as
#' a Horn hypergeometric function or Humbert's hypergeometric used in Gordy
#' (1998) with integral representation:
#'
#'  phi_1(a,b,c,x,y) = Beta(a,b) Int_0^1
#' t^(a-1) (1 - t)^(c-a-1) (1 - yt)^(-b) exp(x t) dt
#' \url{https://en.wikipedia.org/wiki/Humbert_series} Note that Gordy's
#' arguments for x and y are reversed in the reference above.
#'
#' Code for phi1 provided by Gordy.
#'
#'
#' @param a a > 0
#' @param b arbitrary
#' @param c c > 0
#' @param x x > 0
#' @param y 0 <= y < 1
#' @author Merlise Clyde (\email{clyde@@stat.duke.edu})
#' @references Gordy 1998
#' @keywords math
#' @examples
#'
#' # special cases
#' # phi1(a, b, c, x=0, y) is the same as 2F1(b, a; c, y)
#' phi1(1, 2, 1.5, 0, 1 / 100)
#' hypergeometric2F1(2, 1, 1.5, 1 / 100, log = FALSE)
#'
#' phi1(1, 0, 1.5, 3, 1 / 100)
#' hypergeometric1F1(1, 1.5, 3, log = FALSE)
#' @rdname phi1
#' @family special functions
#' @export
#'
#'
phi1 <- function(a, b, c, x, y) {
  # phi1 = int u^{t-1} (1 - v u)^{q - 1} e^-{s u} /B(t, q) (theta
  n <- length(t)
  out <- rep(0, n)
  ans <- .C(C_phi1, as.numeric(a), as.numeric(b), as.numeric(c), as.numeric(x), as.numeric(y), out = as.numeric(out), as.integer(n))$out
  return(ans)
}
