test_that("truncated priors normalization", {
  # Poisson
  # a function for log(1-exp(-x)) 
  # based on Martin M\"{a}chler's work (https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf)
  log1mexp = function(x, x0=log(2)) ifelse(x<x0, log(-expm1(-x)), log1p(-exp(-x)))
  
  
  # generate data
  p = 3
  n = 2*p
  X = matrix(rnorm(n*p),n,p)
  y = rnorm(n)
  
  ### truncated poisson
  # poisson rate
  lambda = 1
  trunc = p-1
  
  # get log prior
  log_prior = dpois(0:trunc, lambda = lambda, log=TRUE) - ppois(trunc, lambda = lambda, log.p = TRUE) - lchoose(p,0:trunc)
  
  # get BAS output
  bas_out = bas.lm(y~X, modelprior = tr.poisson(lambda = lambda, trunc = p-1))
  
  # check prior, something is wrong
  expect_equal(1, sum(bas_out$priorprobs))
  expect_equal(1, sum(exp(log_prior[bas_out$size])))
  expect_equal(bas_out$priorprobs, exp(log_prior[bas_out$size]))
  
  ### truncated power prior

  # power
  kappa = 2
  
  # get log prior
  # uses (p+1)^(-kappa*|gamma|) instead of p^(-kappa*|gamma|) 
  log_prior = (-kappa*c(0:p))*log(trunc +1) - (log1mexp((trunc+1)*kappa*log(trunc+1)) - log1mexp(kappa*log(trunc+1))) - 
               lchoose(p, 0:p)
  
  #get BAS output
  bas_out = bas.lm(y~X, modelprior = tr.power.prior(kappa = kappa, trunc = p-1))
  
  #check prior
  expect_equal(1,sum(bas_out$priorprobs))
  expect_equal(1, sum(exp(log_prior[bas_out$size])))
  expect_equal(bas_out$priorprobs, exp(log_prior[bas_out$size]))
  
})
