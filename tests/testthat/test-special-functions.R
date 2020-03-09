context("special functions")

test_that("phi1", {
  expect_equal(
    phi1(1, 2, 1.5, 0, 1 / 100),
    hypergeometric2F1(2, 1, 1.5, 1 / 100, log = FALSE)
  )
  expect_equal(
    phi1(1, 0, 1.5, 3, 1 / 100),
    hypergeometric1F1(1, 1.5, 3, log = FALSE)
  )
  expect_error(phi1(c(1, 1), c(2, 2), c(1.5, 1.5), c(3, 3), c(.1)))
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
expect_equal(TRUE, !is.na(hypergeometric2F1(12, 1, 2, .85, method = "Laplace")))
})

