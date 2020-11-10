test_that("eplogp-nonfullrank", {
  skip_on_os("solaris")
  loc <- system.file("testdata", package = "BAS")
  d <- read.csv(paste(loc, "JASP-testdata.csv", sep = "/"))

  fullModelFormula <- as.formula("contNormal ~  contGamma * contExpon + contGamma * contcor1 + contExpon * contcor1")


  expect_error(eplogprob(lm(fullModelFormula, data = d)))

})

# issue #54  eplogprob returns NaN

test_that("eplogp-highlysignif", {
  skip_on_os("solaris")
  loc <- system.file("testdata", package = "BAS")
  load(paste0(loc, "/eplogp-testdata.Rdata"))
  prob = eplogprob(lm(Y ~ ., data = df))
  expect_true(sum(is.na(prob)) == 0)

})
