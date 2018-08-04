context("bas.lm")

coverage(codecov)
use_coverage(pkg='.', type=c("codecov"))
test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
