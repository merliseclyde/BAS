context("plot.bas")

test_that("test basic BAS plots", {
  data(Hald)
  hald.gprior <- bas.lm(Y ~ ., data = Hald, prior = "g-prior",
                        modelprior = beta.binomial(1, 1),
                        initprobs = "eplogp")
  expect_null(plot(hald.gprior, drop.always.included=TRUE, ask = FALSE))
  
  expect_error(plot(hald.gprior, which = 5))
  hald.gprior <- bas.lm(Y ~ ., include.always = Y ~ X1 + X4, data = Hald, prior = "g-prior",
                        modelprior = beta.binomial(1, 1),
                        initprobs = "eplogp")
 
  expect_null(plot(hald.gprior, drop.always.included=TRUE, ask=FALSE))
 
})
