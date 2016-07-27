list2matrix.bas = function(x, what, which.models=NULL) {
  namesx = x$namesx
  if (is.null(which.models)) which.models= 1:x$n.models
  
  listobj = x[[what]][which.models]
  which = x$which[which.models]
  n.models = length(which.models)
  p = length(namesx)
  mat = matrix(0, nrow=n.models, ncol=p)

  for (i in 1:n.models) {
    mat[i, which[[i]]+1] = listobj[[i]]
  }
  colnames(mat) = namesx
  return(mat)
}

list2matrix.which = function(x, which.models=NULL) {

    namesx = x$namesx
    listobj = x$which
    if (!is.null(which.models)) listobj = listobj[which.models]
    p = length(namesx)
    mat = t(sapply(listobj,
      function(x, dimp) {
        xx = rep(0,dimp)
        xx[x+1] = 1
        xx},
      p))
    colnames(mat) = namesx
    mat}

 

  which.matrix = function(which, n.vars) {
  mat = t(sapply(which,
                  function(x, dimp) {
                    xx = rep(0,dimp)
                    xx[x+1] = 1
                    xx},
                  n.vars))
  mat}


bin2int = function(model) { 
  if (length(model) > 1) { 
    i = sum(2^(model[-1] - 1)) + 1}
  else{ i = 1}
#  if (!is.integer(i)) warning("Exceeded the largest integer for this machine")
  return(unlist(i))
}




