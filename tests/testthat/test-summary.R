context("summary")
test_that("summary.bas", {
  data(Hald)
  hald.bas <- bas.lm(Y ~ .,
                     prior = "BIC",
                     modelprior = uniform(), data = Hald
  )
  expect_length(summary(hald.bas), 60)
  expect_null(print(hald.bas))
})
