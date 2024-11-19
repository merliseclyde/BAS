skip_on_cran()
test_that("count.heredity.models", {
  
  # for checking correctness, generates all models and then omits those that
  # fail the hereditary check
  no_models_f <- function(f) {
    
    t           <- terms(f)
    all_factors <- attr(t, "factors")[-1, , drop = FALSE]
    o           <- attr(t, "order")
    
    # adjacency matrix such that adj[i, j] == 1 implies that i is a parent of j.
    adj_mat <- crossprod(all_factors != 0)
    for (i in seq_len(nrow(adj_mat))) {
      adj_mat[i, ] <- adj_mat[i, ] == o
    }
    tb_o <- tabulate(o) # how often there is an nth-order interaction
    
    # list of parent indices for each term
    required_lower_order <- lapply(seq_len(ncol(all_factors)), function(i) {
      return(which(adj_mat[i, seq_len(i - 1)] == 1))
    })
    
    # generate all possible models
    all_models <- as.matrix(expand.grid(rep(list(0:1), ncol(all_factors))))
    colnames(all_models) <- colnames(all_factors)
    colIndices <- which(lengths(required_lower_order) > 0)
    
    # loop over all models and remove those that fail the hereditary check
    keep <- rep(TRUE, nrow(all_models))
    for (i in colIndices) {
      
      n <- length(required_lower_order[[i]])
      keep <- keep & (
        (rowSums(all_models[ , required_lower_order[[i]], drop = FALSE]) == n) | (all_models[ , i] == 0)
      )
      
    }
    
    all_models <- all_models[keep, , drop = FALSE]
    all_models
  }
  
  brute_force <- function(f) {
    nrow(no_models_f(f))
  }
  
  # analytic version of brute_force
  analytic <- function(f) {
    
    t           <- terms(f)
    all_factors <- attr(t, "factors")[-1, , drop = FALSE]
    o           <- attr(t, "order")
    
    # base case 1: all main effects
    if (all(o == 1))
      return(2^length(o))
    else {
      
      adj_mat <- crossprod(all_factors != 0)
      for (i in seq_len(nrow(adj_mat))) {
        adj_mat[i, ] <- adj_mat[i, ] == o
      }
      tb_o <- tabulate(o)
      
      # m contains the number of models for each order
      m    <- 0 * c(tb_o)
      m[1] <- 2^tb_o[1] # no. models with only main effects
      
      # loop upward over higher order interactions
      for (i in 2:length(tb_o)) {
        
        no_lower_order <- sum(tb_o[seq_len(i - 1)])
        
        for (j in seq_len(tb_o[i])) {
          
          # enumerate all combinations of interactions
          combs <- utils::combn(tb_o[i], j)
          for (k in seq_len(ncol(combs))) {
            
            # the subset of the adjacency matrix that corresponds to the current combination
            subset <- adj_mat[no_lower_order + combs[, k], 1:no_lower_order, drop = FALSE]
            
            # all edges imply a parent which must be included (i.e., a constraint)
            constraints <- .colSums(subset, m = j, n = no_lower_order) > 0
            
            # if there are missing constraints, we need to check if they involve
            # interactions of various levels. If they do we need to recurse
            missing <- which(constraints == 0)
            if (any(o[missing] > 1) && length(unique(o[missing])) > 1) {
              
              # recursive case : construct the formula corresponding to the constraint and
              # call analytic again
              # this can probably be sped up by only passing (subsets) of the adjacency matrix around
              cnms <- colnames(adj_mat)[1:no_lower_order]
              new_f_str <- paste("y ~", paste(cnms[missing], collapse = " + "))
              # print(paste(deparse(f), "->", new_f_str))
              new_f <- as.formula(paste("y ~", paste(cnms[missing], collapse = " + ")))
              no_models <- Recall(new_f)
              
            } else {
              
              # base case 2: all free variables are effects of the same order (e.g,. all main effects, all 2-way interactions, etc.)
              no_constraints <- sum(constraints)
              no_models <- 2^(no_lower_order - no_constraints)
              
            }
            
            m[i] <- m[i] + no_models
            
            # print(subset)
            # cat(sprintf(
            #   "order = %d, no_combs = %d, comb = %d, no_constraints = %d, no_models = %d\n\n", i, j, k, no_constraints, 2^(no_lower_order - no_constraints)
            # ))
          }
        }
      }
      return(sum(m))
    }
  }

  # see BayesFactor:::enumerateGeneralModels for a way to enumerate all possible models
  # accounting for the restrictions.
  # opts <- BayesFactor:::enumerateGeneralModels(y ~ x1 * x2 * x3 * x4, "withmain")
  # for (f in opts)
  #   cat(deparse(f, 300), ",\n", sep = "")
  fs <- list(
    y ~ x1,
    y ~ x2,
    y ~ x1 + x2,
    y ~ x1 + x2 + x1:x2,
    y ~ x3,
    y ~ x1 + x3,
    y ~ x2 + x3,
    y ~ x1 + x2 + x3,
    y ~ x1 + x2 + x1:x2 + x3,
    y ~ x1 + x3 + x1:x3,
    y ~ x1 + x2 + x3 + x1:x3,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3,
    y ~ x2 + x3 + x2:x3,
    y ~ x1 + x2 + x3 + x2:x3,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3,
    y ~ x4,
    y ~ x1 + x4,
    y ~ x2 + x4,
    y ~ x1 + x2 + x4,
    y ~ x1 + x2 + x1:x2 + x4,
    y ~ x3 + x4,
    y ~ x1 + x3 + x4,
    y ~ x2 + x3 + x4,
    y ~ x1 + x2 + x3 + x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4,
    y ~ x1 + x3 + x1:x3 + x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4,
    y ~ x2 + x3 + x2:x3 + x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4,
    y ~ x1 + x4 + x1:x4,
    y ~ x1 + x2 + x4 + x1:x4,
    y ~ x1 + x2 + x1:x2 + x4 + x1:x4,
    y ~ x1 + x3 + x4 + x1:x4,
    y ~ x1 + x2 + x3 + x4 + x1:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4 + x1:x4,
    y ~ x1 + x3 + x1:x3 + x4 + x1:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4 + x1:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x1:x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4 + x1:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x1:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4,
    y ~ x2 + x4 + x2:x4,
    y ~ x1 + x2 + x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x4 + x2:x4,
    y ~ x2 + x3 + x4 + x2:x4,
    y ~ x1 + x2 + x3 + x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4 + x2:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x2:x4,
    y ~ x2 + x3 + x2:x3 + x4 + x2:x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x2:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x2:x4,
    y ~ x1 + x2 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x3 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4,
    y ~ x1 + x2 + x1:x2 + x4 + x1:x4 + x2:x4 + x1:x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4,
    y ~ x3 + x4 + x3:x4,
    y ~ x1 + x3 + x4 + x3:x4,
    y ~ x2 + x3 + x4 + x3:x4,
    y ~ x1 + x2 + x3 + x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4 + x3:x4,
    y ~ x1 + x3 + x1:x3 + x4 + x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x3:x4,
    y ~ x2 + x3 + x2:x3 + x4 + x3:x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x3:x4,
    y ~ x1 + x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x3 + x1:x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x3:x4,
    y ~ x2 + x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x2:x4 + x3:x4,
    y ~ x2 + x3 + x2:x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x4 + x1:x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4 + x1:x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4 + x1:x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x1:x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4,
    y ~ x1 + x3 + x1:x3 + x4 + x1:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4 + x1:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x1:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4 + x1:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4 + x1:x3:x4,
    y ~ x2 + x3 + x2:x3 + x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x1:x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x1:x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x3:x4 + x1:x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4 + x1:x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4 + x1:x3:x4 + x2:x3:x4,
    y ~ x1 + x2 + x1:x2 + x3 + x1:x3 + x2:x3 + x1:x2:x3 + x4 + x1:x4 + x2:x4 + x1:x2:x4 + x3:x4 + x1:x3:x4 + x2:x3:x4 + x1:x2:x3:x4
  )
  
  # validate that analytic result matches brute force
  expect_equal(sapply(fs, brute_force), sapply(fs, count.heredity.models))
 
  f <- as.formula((paste("y ~", 
                  paste0("x", 1:4, collapse = " * "), " + ", paste0("x", 5:30, collapse = " + "))))
  
  expect_equal(11207180288, count.heredity.models(f))
  
  f <- y ~ x1 + x2 + x3 + x1:x2:x3:x4
  # main effects are present for the 4-way interaction, but
  # no 2-way or 3-way interactions are added
  expect_equal(brute_force(f), count.heredity.models(f))
 
})
