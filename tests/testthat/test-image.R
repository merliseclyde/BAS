context("image.bas")

test_that("test image plots", {
  data("Hald")
  hald.ZSprior <- bas.lm(Y ~ ., data = Hald, prior = "ZS-null")
  expect_null(image(hald.ZSprior, drop.always.included = TRUE))
  expect_null(image(hald.ZSprior, drop.always.included = FALSE))
  expect_null(image(hald.ZSprior, intensity = FALSE))
  expect_null(image(hald.ZSprior, intensity = TRUE, color = "blackandwhite"))
  # tests related to issue #43
  expect_null(image(hald.ZSprior, rotate = FALSE, mar = c(6, 8, 6, 2) + .1))
  expect_null(image(hald.ZSprior, prob = FALSE, rotate = TRUE,
                    mar = c(6, 8, 6, 2) + .1))
  expect_error(image(hald.ZSprior, prob = FALSE, rotate = TRUE,
                     namesx = c("Intercept")))
})
