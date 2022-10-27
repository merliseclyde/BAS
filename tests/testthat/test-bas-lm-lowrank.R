# the following does not throw a warning on a m1MAC  github issue #62 skipping 
# test pending package archival until able to debug on M1mac

test_that("check non-full rank", {
  loc <- system.file("testdata", package = "BAS")
  d <- read.csv(paste(loc, "JASP-testdata.csv", sep = "/"))
  
  fullModelFormula <- as.formula("contNormal ~  contGamma * contExpon + contGamma * contcor1 + contExpon * contcor1")
  
  
  expect_error(eplogprob(lm(fullModelFormula, data = d)))
  
  basObj.eplogp <- bas.lm(fullModelFormula,
                          data = d,
                          alpha = 0.125316, initprobs = "marg-eplogp",
                          prior = "JZS", method = "deterministic", pivot = TRUE,
                          modelprior = uniform(),
                          weights = facFifty, force.heredity = FALSE
  )
  basObj.det <- bas.lm(fullModelFormula,
                       data = d,
                       alpha = 0.125316,
                       modelprior = uniform(),
                       prior = "JZS", method = "deterministic", pivot = TRUE,
                       weights = facFifty, force.heredity = FALSE
  )
  basObj <- bas.lm(fullModelFormula,
                   data = d,
                   alpha = 0.125316, modelprior = uniform(),
                   prior = "JZS", method = "BAS", pivot = TRUE,
                   weights = facFifty, force.heredity = FALSE
  )
  expect_equal(0, sum(is.na(basObj.det$postprobs)))
  expect_equal(basObj.eplogp$probne0, basObj.det$probne0)
  expect_equal(basObj.det$probne0, basObj$probne0)
  expect_equal(basObj.eplogp$probne0, basObj$probne0)
  
  basObj.EBL <- bas.lm(fullModelFormula,
                       data = d,
                       alpha = 0.125316, initprobs = "marg-eplogp",
                       prior = "EB-local", method = "deterministic", pivot = TRUE,
                       modelprior = uniform(),
                       weights = facFifty, force.heredity = FALSE
  )
  basObj.up <- update(basObj.eplogp, newprior = "EB-local")
  expect_equal(basObj.EBL$postprobs, basObj.up$postprobs,
               tolerance=.001)
  
  skip_on_cran()
  expect_warning(bas.lm(fullModelFormula,
                        data = d,
                        alpha = 0.125316,
                        prior = "JZS",
                        weights = facFifty, force.heredity = FALSE, pivot = FALSE))
})
