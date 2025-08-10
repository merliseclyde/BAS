context(" issue 96 rank deficient models and R2")
# Issue #96
test_that("rank and R2", {
  # generate data and get bas output 
  n = 4
  X = cbind(matrix(1, n, n), diag(1, n, n), matrix(1, n, n) - diag(1, n, n))
  y = rnorm(n)
  bas_out = bas.lm(y~X, data=data.frame(y, X))
  cond = bas_out$R2[bas_out$rank < n & n <= bas_out$size] < 1 # should be > 1
  expect_equal(sum(cond), length(cond))
  cond = bas_out$R2[bas_out$rank == n & n <= bas_out$size] == 1
  expect_equal(sum(cond), length(cond))
  })

