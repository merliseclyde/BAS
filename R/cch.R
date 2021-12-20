#' Compound Confluent hypergeometric function of two variables
#'
#' Compute the Confluent Hypergeometric function of two variables, also know as
#' a Horn hypergeometric function or Humbert's hypergeometric used in Gordy
#' (1998) with integral representation:
#'
#'  phi_1(a,b,c,x,y) =  [(Gamma(c)/Gamma(a) Gamma(a-c))] Int_0^1
#' t^(a-1) (1 - t)^(c-a-1) (1 - yt)^(-b) exp(x t) dt
#' \url{https://en.wikipedia.org/wiki/Humbert_series} Note that Gordy's
#' arguments for x and y are reversed in the reference above.
#'
#' The original `phi1` function in `BAS` was based on `C` code provided by 
#' Gordy.  This function returns NA's 
#' when x is greater than `log(.Machine$double.xmax)/2`.   A more 
#' stable method for calculating the `phi1` function using R's `integrate` 
#' was suggested by Daniel Heemann and is now an option whenever $x$ is too 
#' large.  For calculating Bayes factors that use the `phi1` function we 
#' recommend using the `log=TRUE` option to compute log Bayes factors.
#'
#'
#' @param a a > 0
#' @param b arbitrary
#' @param c c > 0
#' @param x x > 0
#' @param y y > 0 
#' @param log logical indicating whether to return phi1 on the log scale
#' @author Merlise Clyde (\email{clyde@@duke.edu})
#' @author Daniel Heemann (\email{df.heemann@@gmail.com})
#' @references Gordy 1998
#' @keywords math
#' @examples
#'
#' # special cases
#' # phi1(a, b, c, x=0, y) is the same as 2F1(b, a; c, y)
#' phi1(1, 2, 1.5, 0, 1 / 100, log=FALSE)
#' hypergeometric2F1(2, 1, 1.5, 1 / 100, log = FALSE)
#'
#' # phi1(a,0,c,x,y) is the same as 1F1(a,c,x)
#' phi1(1, 0, 1.5, 3, 1 / 100)
#' hypergeometric1F1(1, 1.5, 3, log = FALSE)
#' 
#' # use direct integration
#' phi1(1, 2, 1.5, 1000, 0, log=TRUE)
#' @rdname phi1
#' @family special functions
#' @export
#'
#'
phi1 <- function(a, b, c, x, y, log=FALSE) {
# phi_1(a,b,c,x,y) = 
#     Int_0^1 t^(a-1) (1 - t)^(c-a-1) (1 - y t)^(-b) exp(x t) dt/Beta(a, c-a)
  na <- length(a)
  nb <- length(b)
  nc <- length(c)
  nx <- length(x)
  ny <- length(y)
  
# if (any(y < 0 | y >= 1) )stop("y is outside of [0, 1)")
# if (any(x < 0)) stop("x must be >= 0")

  ns = c(na,nb, nc, nx, ny)
  n = max(ns)

  if ((n > 1) && (mean(ns) != n)) {
    stop("length of inputs are not the same")
  }
  
  MV = log(.Machine$double.xmax)/2
  MX = max(x)
  div = ceiling(MX/MV)
  scale = 1/exp(max(0, (MX - MV)/div))

  
  out <- rep(0, n)
  ans <- .C(C_phi1,
            as.numeric(a),
            as.numeric(b),
            as.numeric(c),
            as.numeric(x),
            as.numeric(y),
            as.integer(div),
            as.numeric(scale),
            out = as.numeric(out), as.integer(n))$out
  if (!log) ans=exp(ans)
  return(ans)
}

#' Truncated Compound Confluent Hypergeometric function 
#'
#' Compute the Truncated Confluent Hypergeometric function from Li and Clyde
#'  (2018) which is the normalizing constant in the tcch density of Gordy
#' (1998) with integral representation:
#'
#' tr.cch(a,b,r,s,v,k) =  Int_0^1/v
#'                     u^(a-1) (1 - vu)^(b -1) (k + (1 - k)vu)^(-r) exp(-s u) du
#' 
#' This uses a more 
#' stable method for calculating the normalizing constant using R's `integrate` 
#' function rather than the version in Gordy 1998. For calculating Bayes factors 
#' that use the `trCCH` function we 
#' recommend using the `log=TRUE` option to compute log Bayes factors.
#'
#'
#' @param a a > 0
#' @param b b > 0
#' @param r r  >= 0
#' @param s arbitrary
#' @param v 0 < v 
#' @param k arbitrary
#' @param log logical indicating whether to return values on the log scale; 
#' useful for Bayes Factor calculations
#' @author Merlise Clyde (\email{clyde@@duke.edu})
#' @references Gordy 1998 Li & Clyde 2018
#' @keywords math
#' @aliases trunc.CCH
#' @examples
#'
#' # special cases
#' # trCCH(a, b, r, s=0, v = 1, k) is the same as
#' # 2F1(a, r, a + b, 1 - 1/k)*beta(a, b)/k^r
#' 
#' k = 10; a = 1.5; b = 2; r = 2;  
#' trCCH(a, b, r, s=0, v = 1, k=k) *k^r/beta(a,b)
#' hypergeometric2F1(a, r, a + b, 1 - 1/k, log = FALSE)
#'
#' # trCCH(a,b,0,s,1,1) is the same as 
#' # beta(a, b) 1F1(a, a + b, −s, log=FALSE)
#' s = 3; r = 0; v = 1; k = 1
#' beta(a, b)*hypergeometric1F1(a, a+b, -s, log = FALSE)
#' trCCH(a, b, r, s, v, k)
#' 
#' # Equivalence with the Phi1 function 
#' a = 1.5; b = 3; k = 1.25; s = 400;  r = 2;  v = 1; 
#' 
#' phi1(a, r,  a + b, -s, 1 - 1/k,  log=FALSE)*(k^-r)*gamma(a)*gamma(b)/gamma(a+b)
#' trCCH(a,b,r,s,v,k)
#' @rdname trCCH
#' @family special functions
#' @export
trCCH <- function(a, b, r, s, v, k, log=FALSE) {
  #  phi_1(a,b,c,x,y) =  (Gamma(c)/Gamma(a) Gamma(a-c)) Int_0^1
  # t^(a-1) (1 - t)^(c-a-1) (1 - y t)^(-b) exp(x t) dt
  na <- length(a)
  nb <- length(b)
  nr <- length(r)
  ns <- length(s)
  nv <- length(v)
  nk <- length(k)
  
 # if (any(v <= 0 | v > 1) )stop("v is outside of (0, 1]")
 # if (any(s < 0)) stop("s must be >= 0")
  
  ns = c(na,nb, nr, ns, nv, nk)
  n = max(ns)
  
  if ((n > 1) && (mean(ns) != n)) {
    stop("length of inputs are not the same")
  }
  
  div = 1; scale = 1;
#  MV = log(.Machine$double.xmax)/2
#  MX = max(s)
#  div = ceiling(MX/MV)
#  scale = 1/exp(max(0, (MX - MV)/div))
  
  
  out <- rep(0, n)
  ans <- .C(C_tcch,
            as.numeric(a),
            as.numeric(b),
            as.numeric(r),
            as.numeric(s),
            as.numeric(v),
            as.numeric(k),
            as.integer(div),
            as.numeric(scale),
            out = as.numeric(out), as.integer(n))$out
  if (!log) ans=exp(ans)
  return(ans)
}
