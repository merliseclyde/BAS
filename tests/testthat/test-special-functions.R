context("special functions")

test_that("phi1", {
  expect_equal(
    phi1(1, 2, 1.5, 0, 1 / 100),
    hypergeometric2F1(2, 1, 1.5, 1 / 100, log = FALSE)
  )
  expect_equal(
    phi1(1, 0, 1.5, 3, 1 / 100, log=FALSE),
    hypergeometric1F1(1, 1.5, 3, log = FALSE)
  )
  expect_error(phi1(c(1, 1), c(2, 2), c(1.5, 1.5), c(3, 3), c(.1)))
  expect_equal(TRUE, is.finite(phi1(1, 2, 1.5, 1000, 1/100, log=TRUE)))  # Issue #55
  expect_equal(FALSE, is.finite(phi1(1, 2, 1.5, 1000, 1/100, log=FALSE))) 
  expect_length(
    phi1(c(1, 1), c(2, 2), c(1.5, 1.5), c(3, 3), c(.1, .1)),
    2
  )
})

test_that("2F1", {
expect_warning(hypergeometric2F1(1,0,-1, .5))
expect_warning(hypergeometric2F1(1,1,.5, .5))
expect_warning(hypergeometric2F1(1,1,1.5, 1.5))
expect_warning(hypergeometric2F1(1,1,1.5, -1.5))
expect_warning(hypergeometric2F1(10000,1,.5, .99995))
expect_equal(Inf, hypergeometric2F1(1,1,2, 1.0))
expect_equal(TRUE, hypergeometric2F1(1,1,5, 1.0,
                                    method="Laplace",
                                    log=FALSE)>0)
expect_warning(hypergeometric2F1(3,1,1000, .999))
})



#' 

test_that("tcch", {
  k = 10;  a = 2; b = 3; r = 2;  
  expect_equal(
   # trCCH(a, b, r, s=0, v = 1, k) is the same as
   # 2F1(a, r, a + b, 1 - 1/k)*beta(a, b)/k^r
    trCCH(a, b, r, s=0, v = 1, k=k) *k^r/beta(a,b),
    hypergeometric2F1(a, r, a + b, 1 - 1/k, log = FALSE)
  )
  expect_equal(
    # trCCH(a, b, r, s=0, v = 1, k) is the same as
    # 2F1(a, r, a + b, 1 - 1/k)*beta(a, b)/k^r
    trCCH(a, b, r, s=0, v = 1, k=k, log=TRUE) +log(k)*r - lbeta(a,b),
    hypergeometric2F1(a, r, a + b, 1 - 1/k, log = TRUE)
  )
  s = 3; r = 0; v = 1; k = 10
  # beta(a, b)*hypergeometric1F1(a, a+b, -s, log = FALSE) is the same as 
  # trCCH(a, b, r, s, v, k)
  expect_equal(
    beta(a, b)*hypergeometric1F1(a, a+b, -s, log = FALSE) ,
    trCCH(a, b, r, s, v, k)
  )
  a = 1.5; b = 3; k = 1.25; s = 40;  r = 2;  v = 1; k = 1.25
  expect_equal(
    phi1(a, r,  a + b, -s, 1 - 1/k,  log=FALSE)*(k^-r)*gamma(a)*gamma(b)/gamma(a+b),
    trCCH(a,b,r,s,v,k), tolerance = .00001
  )
  expect_error(trCCH(c(1, 1), c(2, 2), c(1.5, 1.5), c(1, 1), c(1), c(1)))
  s = 10000
  expect_equal(TRUE, is.finite(trCCH(a,b,r,s,v,k, log=TRUE)))  # Issue #55
  expect_equal(TRUE, is.finite(trCCH(a,b,r,s,v,k, log=FALSE))) 
  expect_length(trCCH(c(1, 1), c(2, 2), c(1.5, 1.5), c(1, 1), c(1,1), c(1,1)),
    2
  )
})
