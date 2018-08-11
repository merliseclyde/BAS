context("outliers")

test_that("outliers.R", {
   data("stackloss")
   stack_lm <- lm(stack.loss ~ ., data = stackloss)
   stack_outliers <- Bayes.outlier(stack_lm, k = 3)
   expect_null(plot(stack_outliers$prob.outlier, type = "h",
                    ylab = "Posterior Probability"))
   # adjust for sample size for calculating prior prob that a
   # a case is an outlier
   stack_outliers <- Bayes.outlier(stack_lm, prior.prob = 0.95)
  #' # cases where posterior probability exceeds prior probability
  expect_length(which(stack_outliers$prob.outlier >
                        stack_outliers$prior.prob),2)
})
