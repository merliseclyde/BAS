#' Coerce a BAS list object into a matrix.
#' 
#' Models, coefficients, and standard errors in objects of class 'bas' are
#' represented as a list of lists to reduce storage by omitting the zero
#' entries.  These functions coerce the list object to a matrix and fill in the
#' zeros to facilitate other computations.
#' 
#' \code{list2matrix.bas(x, which)} is equivalent to
#' \code{list2matrix.which(x)}, however, the latter uses sapply rather than a
#' loop.
#' \code{list2matrix.which} and \code{which.matrix} both coerce
#' \code{x$which} into a matrix.
#' 
#' @aliases list2matrix.bas list2matrix.which which.matrix
#' @param x a 'bas' object
#' @param what name of bas list to coerce
#' @param which.models a vector of indices use to extract a subset
#' @param which \code{x$which} a list of lists of model indicators
#' @param n.vars the total number of predictors, \code{x$n.vars}
#' @return a matrix representation of \code{x$what}, with number of rows equal
#' to the length of which.models or total number of models and number of
#' columns \code{x$n.vars}
#' @author Merlise Clyde \email{clyde@@stat.duke.edu}
#' @seealso \code{\link{bas}}
#' @keywords regression
#' 
#' 
#' @examples
#' 
#' \dontrun{library(MASS)
#' data(UScrime)
#' UScrime[,-2] = log(UScrime[,-2])
#' crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",
#'                     initprobs= "eplogp") 
#' coef = list2matrix.bas(crime.bic, "ols")  # extract all ols coefficients
#' se = list2matrix.bas(crime.bic, "ols.se")
#' models = list2matrix.which(crime.bic)     #matrix of model indicators
#' models = which.matrix(crime.bic$which, crime.bic$n.vars)     #matrix of model indicators}
#' 
#' @rdname list2matrix
#' @family as.matrix methods 
#' @export
list2matrix.bas = function(x, what, which.models=NULL) {
  namesx = x$namesx
  if (is.null(which.models)) which.models= 1:x$n.models
  
  listobj = x[[what]][which.models]
  which = x$which[which.models]
  n.models = length(which.models)
  p = length(namesx)
  mat = matrix(0, nrow=n.models, ncol=p)

  for (i in 1:n.models) {
    mat[i, which[[i]]+1] = listobj[[i]]
  }
  colnames(mat) = namesx
  return(mat)
}
#' @rdname list2matrix
#' @family as.matrix methods
#' @export
list2matrix.which = function(x, which.models=NULL) {

    namesx = x$namesx
    listobj = x$which
    if (!is.null(which.models)) listobj = listobj[which.models]
    p = length(namesx)
    mat = t(sapply(listobj,
      function(x, dimp) {
        xx = rep(0,dimp)
        xx[x+1] = 1
        xx},
      p))
    colnames(mat) = namesx
    mat}
#' Coerce a BAS list object of models into a matrix.
#' 
#'  This function coerces the list object of models to a matrix and fill in the
#' zeros to facilitate other computations.
#' 
#' \code{which.matrix}  coerces
#' \code{x$which} into a matrix.
#' 
#' @aliases which.matrix
#' @param which a 'bas' model object  \code{x$which}
#' @param n.vars the total number of predictors, \code{x$n.vars}
#' @return a matrix representation of \code{x$which}, with number of rows equal
#' to the length of which.models or total number of models and number of
#' columns \code{x$n.vars}
#' @author Merlise Clyde \email{clyde@@stat.duke.edu}
#' @seealso \code{\link{bas}}
#' @keywords regression
#' 
#' 
#' @examples
#' 
#' \dontrun{library(MASS)
#' data(UScrime)
#' UScrime[,-2] = log(UScrime[,-2])
#' crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",
#'                     initprobs= "eplogp") 
#' models = which.matrix(crime.bic$which, crime.bic$n.vars)     #matrix of model indicators}
#' 
#' @rdname which.matrix
#' @family as.matrix methods
#' @export
  which.matrix = function(which, n.vars) {
  mat = t(sapply(which,
                  function(x, dimp) {
                    xx = rep(0,dimp)
                    xx[x+1] = 1
                    xx},
                  n.vars))
  mat}



