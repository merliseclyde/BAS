.make.parents.of.interactions =
  function (mf, data) 
  { termnamesX = attr(terms(mf), "term.labels") 
    p = length(termnamesX)
    interactions = grep(":", termnamesX)
    parents = diag(p)
    colnames(parents) = termnamesX
    rownames(parents) = termnamesX
    for (j in interactions) {
      term = interactions[j]
      main = unlist(strsplit(termnamesX[j], ":", 
                                      fixed = TRUE))
        parents.of.term = main
       for (i in 2:length(main)) {
        parents.of.term = c(parents.of.term, 
                            combn(main, i, FUN=paste0, collapse=":"))
              }
      parents[j, parents.of.term] = 1
    }
    
    X = model.matrix(mf, data)
    loc = attr(X, "assign")[-1] #drop intercept
    parents = parents[loc,]
    row.names(parents) = colnames(X)[,-1]
    return(list(X = X, parents=parents))
  }


# model.matrix(mf, data) 
# attr( , "assign") has where terms are located

# mp = .make.parents.of.interactions(mf, df)

