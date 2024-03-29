context("image.bas")

test_that("test image plots", {
  data("Hald")
  hald.ZSprior <- bas.lm(Y ~ ., data = Hald, prior = "JZS")
  expect_null(image(hald.ZSprior, drop.always.included = TRUE))
  expect_null(image(hald.ZSprior, drop.always.included = FALSE))
  expect_null(image(hald.ZSprior, rotate = FALSE))
  expect_null(image(hald.ZSprior, prob = FALSE))
  expect_null(image(hald.ZSprior, intensity = FALSE))
  expect_null(image(hald.ZSprior, intensity = TRUE, color = "blackandwhite"))
  hald.ZSprior <- bas.lm(Y ~ ., data = Hald, include.always = Y ~ ., prior = "JZS")
  expect_error(image(hald.ZSprior, drop.always.included=TRUE))
})
