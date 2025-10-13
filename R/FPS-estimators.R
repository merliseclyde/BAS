# FPS Estimators
# 
compute_sample_probs_Bayes_HT = function(object) {
  models <- which.matrix(object$which, object$n.vars)
  probs <- object$probne0
  exp(models %*% log(probs) + (1 - models)%*% log(1 - probs))
}