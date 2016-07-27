cv.summary.bas = function(pred, ytrue, score="squared-error") {
  if (length(pred) != length(ytrue)) {
    warning("predicted values and observed values are not the same length")
    return()
  }
  if (score == "miss-class") {
    pred.class <- ifelse(pred < .5, 0, 1)
    confusion = table(pred.class, ytrue)
    error = (sum(confusion) - sum(diag(confusion)))/sum(confusion) 
  }
  else  {
     error = sqrt(sum((pred - ytrue)^2)/length(ytrue))
  }
  return(error)
}
