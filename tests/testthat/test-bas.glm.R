test_that("GLM logit probne0", {
  data(Pima.tr, package="MASS")
  pima.BAS = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
                     method="BAS",
                     betaprior=bic.prior(), family=binomial(),
                     modelprior=uniform())
  pima.det = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
                     method="deterministic",
                     betaprior=bic.prior(), family=binomial(),
                     modelprior=uniform())
  expect_equal(pima.BAS$probne0, pima.det$probne0)
  expect_equal(predict(pima.BAS, type='link', se.fit=TRUE)$fit,
               as.vector(fitted(pima.det)))
  expect_equal(predict(pima.BAS, type='response')$fit,
               fitted(pima.det, type='response'))
  pima.BAS = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
                     method="BAS",
                     betaprior=Jeffreys(), family=binomial(),
                     modelprior=uniform())
  pima.det = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
                     method="deterministic",
                     betaprior=Jeffreys(), family=binomial(),
                     modelprior=uniform())
  expect_equal(pima.BAS$probne0, pima.det$probne0)

})
