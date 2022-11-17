context("update.bas")

test_that("update bas.lm", {
  data(Hald)
  hald_bic <- bas.lm(Y ~ .,
                     data = Hald, prior = "BIC",
                     initprobs = "eplogp", method="deterministic")
  
  expect_error(update(hald_bic, newprior="g-prior"))
  
  hald_EBG <- bas.lm(Y ~ .,
                     data = Hald, prior = "EB-global",
                     initprobs = "eplogp", method="deterministic")
  hald_update <- update(hald_bic, newprior="EB-global")
  
  expect_equal(hald_EBG$postprobs, hald_update$postprobs)
}
)