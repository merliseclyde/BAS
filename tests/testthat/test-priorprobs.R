# issue #87 Prior inclusions probabilities appear incorrect for a Bernoulli(.2) 
# prior when always including one other predictor


test_that("bernoulli prior and include.always", {

# helper function to compute prior inclusion probabilities
compute_prior_inclusion_probs <- function(fit) {
 
  modelProbs <- fit$priorprobs/sum(fit$priorprobs)
  priorProbs <- as.vector(t(modelProbs)%*% which.matrix(fit$which, fit$n.vars))
  names(priorProbs) <- fit$namesx
  
  return(priorProbs)
}

# some basic data
n <- 100
p <- 6
x <- matrix(rnorm(n*p), n, p)
beta <- rnorm(p)
y <- x %*% beta + rnorm(n)

data <- data.frame(y, x)
colnames(data) <- c("y", letters[1:p])

# wrong 
fit1 <- bas.lm(
  formula         = y ~ .,
  data            = data,
  prior           = "hyper-g-laplace",
  alpha           = 3,
  modelprior      = Bernoulli(.2),
  n.models        = NULL,
  method          = "BAS",
  MCMC.iterations = NULL,
  include.always = ~ a,
#  initprobs       = c(1, 1, seq(.5, 0.1, length=5)),
  weights         = NULL,
  renormalize     = TRUE
) 



# wrong still
fit2 <- bas.lm(
  formula         = y ~ .,
  data            = data,
  prior           = "hyper-g-laplace",
  alpha           = 3,
  modelprior      = Bernoulli(c(1,1, rep(.2,5))),
  n.models        = NULL,
  method          = "BAS",
  MCMC.iterations = NULL,
  include.always = ~ a,
  #  initprobs       = c(1, 1, seq(.5, 0.1, length=5)),
  weights         = NULL,
  renormalize     = TRUE
) 

# as expected 

fit3 <- bas.lm(
  formula         = y ~ .,
  data            = data,
  prior           = "hyper-g-laplace",
  alpha           = 3,
  modelprior      = Bernoulli(c(1,1, rep(.2,5))),
#  include.always = ~ a,                               
  n.models        = NULL,
  method          = "BAS",
  MCMC.iterations = NULL,
#  initprobs       = c(1, seq(.5, 0.1, length=5)),
  weights         = NULL,
  renormalize     = TRUE
) 

expect_equal(compute_prior_inclusion_probs(fit1), compute_prior_inclusion_probs(fit2))

hyp = c(1,.5, 1, .5, 1, .5, .5)

fit4 <- bas.lm(
  formula         = y ~ .,
  data            = data,
  prior           = "hyper-g-laplace",
  alpha           = 3,
  modelprior      = Bernoulli(hyp),
  n.models        = NULL,
  method          = "BAS",
  MCMC.iterations = NULL,
  #  initprobs       = c(1, seq(.5, 0.1, length=5)),
  weights         = NULL,
  renormalize     = TRUE
) 

expect_equal(hyp, as.numeric(compute_prior_inclusion_probs(fit4)))

# issue #87 should be equal
expect_equal(compute_prior_inclusion_probs(fit2), compute_prior_inclusion_probs(fit3))
}
)


