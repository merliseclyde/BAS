context("diagnostics")

test_that("diagnostics.R", {
 data(Hald)
 hald_gprior <- bas.lm(Y ~ ., data=Hald, method='MCMC', MCMC.iterations = 10^6)
 expect_null(diagnostics(hald_gprior))
})
