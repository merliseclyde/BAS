test_that("test basic BAS plots", {
  data(Hald)
  hald.gprior =  bas.lm(Y~ ., data=Hald, prior="g-prior", alpha=13,
                        modelprior=beta.binomial(1,1),
                        initprobs="eplogp")
  expect_null(plot(hald.gprior, ask=FALSE))

})
