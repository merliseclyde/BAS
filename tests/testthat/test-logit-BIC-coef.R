test_that("coef in logistic glm", {
  Nrep = 1000
  
  # We load dataset
  data("birthwt", package = "MASS")
  
  birthwt$race <- as.factor(birthwt$race)
  
  # bas.glm
  set.seed(1000000)
  birthwt.bas.glm <- bas.glm(low ~ . -bwt, data=birthwt,
                             family=binomial(link = "logit"),
                             method="MCMC",
                             MCMC.iterations=Nrep,
                             laplace=TRUE)
  expect_no_error(coef(birthwt.bas.glm))
})
