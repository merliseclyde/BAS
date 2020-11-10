test_that("extract Median Probability Model", {

  data(Hald, package="BAS")
  hald_bic =  bas.lm(Y ~ ., data=Hald, alpha=13, prior="BIC",
                     modelprior = uniform())
  hald_MPM_manual =  bas.lm(Y ~ ., data=Hald, alpha=13, prior="BIC",
                            modelprior = uniform(),
                            n.models = 1L,
                            bestmodel = as.numeric(hald_bic$probne0 > .5)
                            )
  hald_MPM = extract_MPM(hald_bic)
  expect_equal(hald_bic$n.vars, hald_MPM$n.vars)
  expect_equal(as.numeric(hald_bic$probne0 > .5),
               as.vector(which.matrix(hald_MPM$which[1], hald_MPM$n.vars)))
  expect_equal(predict(hald_bic, estimator="MPM")$fit,
               predict(hald_MPM)$fit,
               check.attributes = FALSE)

  data(Pima.tr, package="MASS")
  Pima_bas = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7, method="BAS",
                     betaprior=CCH(a=1, b=nrow(Pima.tr)/2, s=0),
                     family=binomial(),
                     modelprior=uniform())
  Pima_MPM_man = bas.glm(type ~ ., data=Pima.tr,  method="BAS",
                        betaprior=CCH(a=1, b=nrow(Pima.tr)/2, s=0),
                        family=binomial(),
                        modelprior=uniform(),
                        n.models = 1L,
                        bestmodel = Pima_bas$probne0 > 0.5)


  Pima_MPM = extract_MPM(Pima_bas)
  expect_equal(as.numeric(Pima_bas$probne0 > .5),
               as.vector(which.matrix(Pima_MPM$which[1], Pima_MPM$n.vars)))
  expect_equal(coef(Pima_MPM)$coef, coef(Pima_MPM_man)$coef)
})
