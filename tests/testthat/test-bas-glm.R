context("bas.glm")

test_that("bas.glm initprobs" , {
 data(Pima.tr, package="MASS")
 expect_error(bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      initprobs=rep(.4, nrow(Pima.tr)-1),
                      betaprior = bic.prior(), family = binomial(),
                      modelprior = uniform())
  )
  set.seed(1)
  pima_bas1 <- bas.glm(type ~ .,
                     data = Pima.tr, method = "BAS",
                     initprobs=rep(.4, ncol(Pima.tr)-1),
                     betaprior = bic.prior(), family = binomial(),
                     modelprior = uniform())
  set.seed(1)
  pima_bas2 <- bas.glm(type ~ .,
                     data = Pima.tr, method = "BAS",
                     initprobs=c(1,rep(.4, ncol(Pima.tr)-1)),
                     betaprior = bic.prior(), family = binomial(),
                     modelprior = uniform())
  expect_equal(pima_bas1$postprobs, pima_bas2$postprobs)
  set.seed(1)
  pima_bas2 <- bas.glm(type ~ .,
                       data = Pima.tr, method = "BAS",
                       initprobs=c(1.5,rep(.4, ncol(Pima.tr)-1)),
                       betaprior = bic.prior(), family = binomial(),
                       modelprior = uniform())
  expect_equal(pima_bas1$postprobs, pima_bas2$postprobs)
})

test_that("GLM logit", {
  data(Pima.tr, package = "MASS")
  set.seed(1)
  pima_BAS <- bas.glm(type ~ .,
    data = Pima.tr, method = "BAS",
    betaprior = bic.prior(), family = binomial(),
    modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$probne0[-1] > 1))
  set.seed(1)
  pima_BAS2 <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = bic.prior(), family = "binomial",
                      modelprior = uniform()
  )
  expect_equal(pima_BAS$postprobs, pima_BAS2$postprobs)
  expect_error(bas.glm(type ~ .,
                       data = Pima.tr, method = "BAS",
                       betaprior = bic.prior(),
                       family = "gaussian",
                       modelprior = uniform())
  )
  expect_error(bas.glm(type ~ .,
                       data = Pima.tr, method = "BAS",
                       betaprior = bic.prior(),
                       family = "homeless",
                       modelprior = uniform())
  )
  pima_det <- bas.glm(type ~ ., data = Pima.tr,
    method = "deterministic", betaprior = bic.prior(),
    family = binomial(), modelprior = uniform()
  )
  expect_equal(pima_BAS$probne0, pima_det$probne0)
  expect_equal(
    predict(pima_BAS, type = "link")$fit,
    as.vector(fitted(pima_det))
  )
  expect_equal(
    predict(pima_BAS, data = Pima.tr, type = "link", se.fit = TRUE)$se.bma.fit,
    predict(pima_BAS, type = "link", se.fit = TRUE)$se.bma.fit
  )
  expect_equal(
    predict(pima_BAS, data = Pima.tr, type = "response", se.fit = TRUE)$se.bma.fit,
    predict(pima_BAS, type = "response", se.fit = TRUE)$se.bma.fit
  )
  expect_equal(
    predict(pima_BAS, type = "response")$fit,
    fitted(pima_det, type = "response")
  )
  pima_BAS <- bas.glm(type ~ .,
    data = Pima.tr, method = "BAS",
    betaprior = Jeffreys(), family = binomial(),
    modelprior = tr.beta.binomial(1, 1, 4))
  pima_det <- bas.glm(type ~ ., data = Pima.tr,
    method = "deterministic", betaprior = Jeffreys(),
    family = binomial(), modelprior = tr.beta.binomial(1, 1, 4))
  expect_equal(pima_BAS$probne0, pima_det$probne0)

  pima_BAS <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = Jeffreys(), family = binomial(),
                      modelprior = tr.poisson(2, 4))
  pima_det <- bas.glm(type ~ ., data = Pima.tr,
                      method = "deterministic", betaprior = Jeffreys(),
                      family = binomial(),
                      modelprior = tr.poisson(2, 4))
  expect_equal(pima_BAS$probne0, pima_det$probne0)
})

# issue "Error in bas.lm: $ operator is invalid for atomic vectors" #5

test_that("model prior string", {
  data(Hald)
  expect_error(bas.glm(type ~ ., data = Pima.tr,
                       method = "deterministic", betaprior = bic.prior(),
                       family = binomial(), modelprior = "uniform")
  )
})


test_that("beta prior string", {
  data(Hald)
  expect_error(bas.glm(type ~ ., data = Pima.tr,
                       method = "deterministic", betaprior = "BIC",
                       family = binomial(), modelprior = uniform())
  )
})

test_that("missing data arg", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = bic.prior(), family = binomial(),
                      modelprior = uniform())
  attach(Pima.tr)
  pima_no_data <- bas.glm(type ~ npreg + glu + bp + skin + bmi + ped + age,
                          method = "BAS",
                          betaprior = bic.prior(),
                          family = binomial(),
                          modelprior = uniform())
  expect_equal(pima_BAS$probne0, pima_no_data$probne0)
  })

test_that("poisson regression", {
  data(crabs, package = "glmbb")
  crabs.bas <- bas.glm(satell ~ color * spine * width + weight,
    data = crabs,
    family = poisson(),
    betaprior = EB.local(), modelprior = uniform(),
    method = "MCMC", n.models = 1024, MCMC.iterations = 10000,
    prob.rw = .95
  )
  expect_null(plot(crabs.bas))
  expect_equal(0, sum(crabs.bas$shrinkage > 1))
})

test_that("glm_fit", {
  data(Pima.tr, package = "MASS")
  Y <- as.numeric(Pima.tr$type) - 1
  X <- cbind(1, as.matrix(Pima.tr[, 1:7]))
  pima_new <- bayesglm.fit(X, Y,
    family = binomial(),
    coefprior = bic.prior(n = length(Y))
  )
  pima_orig <- glm(type ~ ., family = binomial(), data = Pima.tr)
  expect_equivalent(pima_new$coefficients, coef(pima_orig))
  expect_equivalent(pima_new$se, summary(pima_orig)$coef[, 2])
  pima_nowts <- bayesglm.fit(X, Y,
    weights = NULL, offset = NULL,
    family = binomial(),
    coefprior = bic.prior(n = length(Y))
  )
  expect_equal(
    pima_new$coefficients,
    pima_nowts$coefficients
  )
})

test_that("robust prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = robust(), family = binomial(),
    modelprior = uniform()
  )
  expect_equal(nrow(Pima.tr), pima_BAS$betaprior$hyper.parameters$n)
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
  
  n = nrow(Pima.tr)
  p = pima_BAS$size - 1
  W = pima_BAS$Q 
  ns = length(W)
  a = rep(1, ns); b  = rep(2, ns);
  r = rep(1.5, ns);  s = rep(0, ns);  v = (n + 1)/(p + 1); k = rep(1, ns)
  shrinkage  = 1 - exp(trCCH((a + p + 2)/2, b/2, r, (s + W)/2, v, k, log=TRUE) - trCCH((a + p)/2, b/2, r, (s + W)/2, v, k, log=TRUE))
  shrinkage[p == 0] = 1
  expect_equal(shrinkage, pima_BAS$shrinkage , tol=.00001)
  
})

test_that("intrinsic prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ .,
    data = Pima.tr, method = "BAS",
    betaprior = intrinsic(n = nrow(Pima.tr)),
    family = binomial(),
    modelprior = uniform()
  )
  pima_BAS_no_n <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = intrinsic(),
                      family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
  expect_equal(0, sum(pima_BAS_no_n$shrinkage > 1))
  expect_equal(pima_BAS$probne0, pima_BAS_no_n$probne0)
})

test_that("TestBF prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = testBF.prior(g = nrow(Pima.tr)),
    family = binomial(),
    modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
})

test_that("hyper.g.n prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
                      betaprior = hyper.g.n(),
                      family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
})

test_that("hyper.g prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
                      betaprior = hyper.g(),
                      family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
  # value of alpha should be greater than 2
  expect_error(bas.glm(type ~ ., data = Pima.tr, method = "BAS",
                       betaprior = hyper.g(alpha=2.0),
                       family = binomial(),
                       modelprior = uniform())
  )
})

# code coverage with Laplace 
test_that("hyper.g prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
                      betaprior = CCH(alpha=3,beta=1), laplace= TRUE,
                      family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
 
})

# FIXED Issue #31
test_that("g/IC prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = g.prior(g = nrow(Pima.tr)),
    family = binomial(),
    modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
                      betaprior = IC.prior(nrow(Pima.tr)),
                      family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
})



test_that("cch prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_cch <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = CCH(2, 2), family = binomial(),
    modelprior = uniform()
  )
  pima_TG <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = TG(), family = binomial(),
    modelprior = uniform()
  )
  expect_equal(pima_cch$probne0, pima_TG$probne0)
})

test_that("TCCHprior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_Tcch <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = tCCH(alpha = 2, b = 2), family = binomial(),
    modelprior = uniform()
  )
  pima_cch <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = CCH(2, 2), family = binomial(),
    modelprior = uniform()
  )

  
  expect_equal(pima_cch$probne0, pima_Tcch$probne0, tolerance = .00001)
  expect_equal(sort(pima_cch$shrinkage), sort(pima_Tcch$shrinkage), tolerance = .001)
  
  pima_tcch <- bas.glm(type ~ ., data = Pima.tr,
                       method = "BAS",
                       betaprior = tCCH(alpha = 1,
                                        beta = 1,
                                        s = .5),
                       family = binomial(),
                       modelprior = uniform()
  )
  expect_equal(0, sum(pima_tcch$shrinkage > 1))
})

test_that("IC.prior", {
  data(Pima.tr, package = "MASS")
  pima_bic <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = bic.prior(), family = binomial(),
                      modelprior = uniform())
  pima_ic <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = IC.prior(log(nrow(Pima.tr))), family = binomial(),
                      modelprior = uniform())
  expect_equal(pima_bic$probne0, pima_ic$probne0)
  })

test_that("cv.summary", {
  data(Pima.tr, package = "MASS")
  data(Pima.te, package = "MASS")
  pima_bic <- bas.glm(type ~ .,
                      data = Pima.tr, method = "MCMC",
                      MCMC.iterations = 10000,
                      betaprior = bic.prior(), family = binomial(),
                      modelprior = uniform())
  pima_pred <- predict(pima_bic, newdata=Pima.te, type="response")
  expect_equal(TRUE, cv.summary.bas(pima_pred$fit,
                               as.numeric(Pima.te$type)) > 0)
  expect_equal(TRUE, cv.summary.bas(pima_pred$fit,
                                 as.numeric(Pima.te$type),
                                 score="miss-class")
               <1)
  expect_error(cv.summary.bas(pima_pred$fit,
                                    as.numeric(Pima.te$type),
                                    score="percent-explained"))
  expect_error(cv.summary.bas(pima_pred$fit,
                              as.numeric(Pima.te$type[-1]),
                              score="percent-explained"))
})

# FIXED issue #28
test_that("diagnostic plot for glm MCMC", {
  data(Pima.tr, package="MASS")
  pima_MCMC <- bas.glm(type ~ .,
                       data = Pima.tr, MCMC.iterations = 1024,
                       method = "MCMC", betaprior = aic.prior(),
                       family = binomial(),
                       modelprior = tr.poisson(2,5))
  expect_null(diagnostics(pima_MCMC, type = "model"))
  expect_null(diagnostics(pima_MCMC, type = "pip"))
})

# FIXED issue #29
test_that("beta.prime prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
                      betaprior = beta.prime(n =nrow(Pima.tr)), family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(as.numeric(nrow(Pima.tr)), pima_BAS$betaprior$hyper.parameters$n)
  pima_BAS_def <- bas.glm(type ~ ., data = Pima.tr,
                          method = "BAS",
                          betaprior = beta.prime(),
                          family = binomial(),
                          modelprior = uniform())
  expect_equal(pima_BAS_def$probne0, pima_BAS$probne0)
  expect_equal(as.numeric(nrow(Pima.tr)),
               pima_BAS_def$betaprior$hyper.parameters$n)
})

# FIXED issue #33
test_that("Jeffreys & MCMC", {
 data(Pima.tr, package="MASS")
 pima_BAS <-  bas.glm(type ~ .,
                      data = Pima.tr, method = "MCMC",
                      betaprior = Jeffreys(),
                      family = binomial(),
                      modelprior = tr.beta.binomial(1, 1, 4))
expect_equal(0, sum(pima_BAS$probne0 > 1))
expect_length(pima_BAS$probne0, ncol(Pima.tr))
})

#  issue #34
test_that("include always MCMC", {
  data("Pima.tr", package="MASS")
  pima_BAS = bas.glm(type ~ .,
                       data = Pima.tr, method = "MCMC",
                       include.always = ~ bp,
                       betaprior = g.prior(g=100), family = binomial(),
                       modelprior = beta.binomial(1, 1))
  x = pima_BAS$probne0[match(c("Intercept", "bp") ,pima_BAS$namesx)]
  expect_equal(2, sum(x), tolerance=.002)

#  pima_BAS = bas.glm(type ~ .,
#                     data = Pima.tr, method = "BAS",
#                     include.always = ~ bp,
#                     betaprior = g.prior(g=100), family = binomial(),
#                     modelprior = beta.binomial(1, 1))
#  expect_equal(2L, sum(pima_BAS$probne0 >= (1.0 - 10*.Machine$double.eps)))
##  check why method='BAS' does not have 1.0 for keep.
})

# FIXED issue #35
test_that("MCMC+BAS: missing MCMC.iterations and n.models arg", {
  data(Pima.tr, package = "MASS")
  set.seed(1)
  pima_BAS <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = bic.prior(),
                      family = binomial(),
                      modelprior = uniform())
  set.seed(1)
  pima_1 <- bas.glm(type ~ ., data=Pima.tr,
                    method = "MCMC+BAS",
                    betaprior = bic.prior(),
                    family = binomial(),
                    modelprior = uniform())
  set.seed(1)
  pima_2 <- bas.glm(type ~ ., data = Pima.tr,
                    method = "MCMC+BAS",
                    betaprior = bic.prior(),
                    n.models=2^7,
                    MCMC.iterations=10000, #default
                    family = binomial(),
                    modelprior = uniform())
  expect_equal(pima_1$probne0, pima_2$probne0)
  expect_equal(pima_BAS$probne0, pima_2$probne0)
  expect_equal(pima_BAS$n.models, pima_1$n.models)
  expect_equal(pima_BAS$n.models, pima_2$n.models)
})

# issue 38 (in progress)  check that it works with other prior
# with prior probabilities; i.e. failed with Jeffreys

test_that("herdity and BAS", {
  skip_on_os("solaris")
  data(Pima.tr, package="MASS")
  pima_BAS <-  bas.glm(type ~ (bp + glu + npreg)^2,
                       data = Pima.tr, method = "BAS",
                       betaprior = bic.prior(),
                       family = binomial(), update=NULL,
                       modelprior =uniform(),
                       force.heredity=TRUE)
  pima_BAS_no <-  bas.glm(type ~ (bp + glu + npreg)^2,
                       data = Pima.tr, method = "BAS",
                       betaprior = bic.prior(),
                       family = binomial(),  update=NULL,
                       modelprior =uniform(),
                       force.heredity=FALSE)
  pima_BAS_no <- force.heredity.bas(pima_BAS_no)
  expect_equal(0L, sum(pima_BAS$probne0 > 1.0))
  expect_equal(0L, sum(pima_BAS_no$probne0[-1] > 1.0))
  expect_equal(pima_BAS$probne0, pima_BAS_no$probne0)
  expect_equal(0L, sum(duplicated(pima_BAS$which)))
})

# issue 55 in progress
test_that("phi1 and NAs in bas.glm", {
  # parameters for the hyper g/n function
  a1 = 1
  b1 = 2
  r1 = 1.5
  s1 = 0
  v1 = 1
  
  example_df_large <- structure(list(Var1 = structure(c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L), 
                                                      .Label = c("A", "B"), class = "factor"), 
                                     Var2 = structure(c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L),
                                                      .Label = c("A", "B"), class = "factor"), 
                                     Var3 = structure(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                                                      .Label = c("A", "B"), class = "factor"), 
                                     Freq = c(120L, 85L, 266L, 301L,  101L, 146L, 523L, 958L)), 
                                class = "data.frame", row.names = c(NA,-8L))
  b <-  bas.glm(Freq ~ Var1 + Var2 + Var3 + Var1*Var3 +Var2*Var3 + Var1*Var2, 
                data=example_df_large,
                family=poisson(),
                betaprior= hyper.g.n(), modelprior=uniform(),
                include.always = "~ 1 + Var1 + Var2 + Var3 + Var1*Var3 + Var2*Var3",
                n.models=2^10, MCMC.iterations=10,
                prob.rw=.95)
  expect_equal(TRUE, is.finite(exp(b$logmarg[2] - b$logmarg[1])))
})


# github issue 61 
 test_that("Jeffreys prior and include.always", {
  data(Pima.tr, package="MASS"); 
  formula <- type ~1 + npreg + glu + bp + bmi + ped; 
  covariates <- ~1 + npreg; 
  # Do not expect error  so second arg is NA
  expect_error(bas.glm(formula = formula, data = Pima.tr, 
                        family = binomial(), 
                        laplace = FALSE, 
                        betaprior = Jeffreys(), 
                        modelprior = uniform(), 
                        method = "BAS", 
                        include.always = covariates ),
               NA)
})

# test BIC priors
test_that("regression coef and IC priors", {
  data(Pima.tr, package="MASS"); 
  formula <- type ~1 + npreg + glu + bp + bmi + ped; 
  covariates <- ~ 1 + npreg + glu + bp + bmi + ped;
  pima.bas = bas.glm(formula = formula, data = Pima.tr, 
                       family = binomial(), 
                       laplace = FALSE, 
                       betaprior = bic.prior(), 
                       modelprior = uniform(), 
                       include.always = covariates,
                       method = "BAS")

  pima.glm =  glm(formula = formula, data = Pima.tr, 
                   family = binomial()) 

  # postmode and MLE under full models should be equal
  expect_equal(as.numeric(coef(pima.bas)$postmean), as.numeric(coef(pima.glm)))
  
  pima.bas = bas.glm(formula = type ~ bp + bmi, data = Pima.tr, 
                     family = binomial(), 
                     betaprior = bic.prior(), 
                     modelprior = uniform(), 
                     method = "BAS")
  
  # github issue
  expect_warning(coef(pima.bas), NA)
  
})

# test gamma model
test_that("gamma regression coef", {
  
  data(wafer, package="faraway")
  wafer_glm <- glm(formula = resist ~ .,
                   family  = Gamma(link = "log"),
                   data    = wafer)
  # postmode and MLE under full models should be equal
  wafer_bas = bas.glm(resist~ ., data=wafer,  include.always = ~ .,
                      betaprior = bic.prior() ,family = Gamma(link = "log"))
  expect_equal(as.numeric(coef(wafer_bas)$postmean), as.numeric(coef(wafer_glm)))
  
  # expect error  but due to glm not bas as link =logit not possible
  # add error checking for BAS
  expect_error(wafer_bas = bas.glm(resist~ ., data=wafer, 
                      betaprior = bic.prior() ,family = Gamma(link = "logit")))
  
  wafer_bas = bas.glm(resist~ ., data=wafer, 
                                     betaprior = bic.prior() ,family = Gamma(link = "log"))
  
  # do not expect warning FIXME
  
  expect_warning(coef(wafer_bas), NA)
  
})
