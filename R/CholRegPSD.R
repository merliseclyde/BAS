CholRegPSD = function(Y, X) {
  XtX = crossprod(X, X)
  cholQ = chol(XtX, pivot=TRUE)
  pivot = attr(cholQ, "pivot")
  r = attr(cholQ, "rank")
  p = ncol(X)
  if (r < p) warning("rank deficient")

  xty = crossprod(X,Y)[pivot]


 alpha = rep(0, p)

#  Solve Q^T alpha = xty
  if (r >= 1)   alpha[1] = xty[1]/cholQ[1,1]
  if (r > 1) {
    for (i in 2:r) {
      alpha[i] = (xty[i]-sum( cholQ[1:(i-1),i]*alpha[1:(i-1)]))/cholQ[i,i]
    }
  }


# solve L beta = alpha

beta = rep(NA, p)
beta[r] = alpha[r]/cholQ[r,r]

if (r > 1) {
  for (i in (r-1):1) {
    beta[i] = (alpha[i] - sum(cholQ[i, (i+1):r]*beta[(i+1):r]))/cholQ[i,i]
  }}

beta = beta[pivot]
beta[is.na(beta)] = 0

#browser()

return(beta)
}
