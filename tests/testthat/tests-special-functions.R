context("special functions")

test_that("phi1", {
expect_equal(phi1(1, 2, 1.5, 0, 1 / 100),
             hypergeometric2F1(2, 1, 1.5, 1 / 100, log = FALSE))
expect_equal(phi1(1, 0, 1.5, 3, 1 / 100),
             hypergeometric1F1(1, 1.5, 3, log = FALSE))

})
