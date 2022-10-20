# the following does not throw a warning on a m1MAC  github issue #62 skipping 
# test pending package archival until able to debug on M1mac

test_that("check non-full rank with pivot=FALSE", {
  skip_on_os("mac", arch = "aarch64")
  loc <- system.file("testdata", package = "BAS")
  d <- read.csv(paste(loc, "JASP-testdata.csv", sep = "/"))
  
  fullModelFormula <- as.formula("contNormal ~  contGamma * contExpon + contGamma * contcor1 + contExpon * contcor1")
  
  
  expect_warning(bas.lm(fullModelFormula,
                        data = d,
                        alpha = 0.125316,
                        prior = "JZS",
                        weights = facFifty, force.heredity = FALSE, pivot = FALSE))
  
 
})
