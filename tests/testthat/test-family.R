context("bas.glm")

# issue #92 
# 
test_that("family links", {
  data(Pima.tr, package="MASS")
  expect_error(bas.glm(type ~ .,
                       data = Pima.tr, method = "BAS",
                       betaprior = bic.prior(), family = binomial(link = "inverse"),
                       modelprior = uniform())
  )
  data(crabs, package = "glmbb")
  expect_error(crabs.bas <- bas.glm(satell ~ color * spine * width + weight,
                                       data = crabs,
                                       family = poisson(link="logit"),
                                       betaprior = EB.local(), modelprior = uniform(),
                                       method = "MCMC", n.models = 1024, MCMC.iterations = 1000
                  ))
  expect_error(crabs.bas <- bas.glm(satell ~ color * spine * width + weight,
                                    data = crabs,
                                    family = poisson(link="identity"),
                                    betaprior = EB.local(), modelprior = uniform(),
                                    method = "MCMC", n.models = 1024, MCMC.iterations = 1000
  ))
  data(wafer, package="faraway")
  expect_error(bas.glm(formula = resist ~ ., include.always = ~ .,
                         betaprior = bic.prior() ,
                         family  = Gamma(link = "identity"),
                         data    = wafer))
  expect_error(bas.glm(formula = resist ~ ., include.always = ~ .,
                         betaprior = bic.prior() ,
                         family  = gaussian(link = "identity"),
                         data    = wafer))
})
