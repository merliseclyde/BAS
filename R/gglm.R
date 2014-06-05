gglm.fit <- function (x, y, 
		weights = rep(1, nobs), 
	    	offset = rep(0, nobs), 
		family = binomial(),
            control = glm.control()) {
  x <- as.matrix(x)
  ynames <- if (is.matrix(y))     rownames(y)
            else names(y)
  nobs <- NROW(y)
  nvars <- ncol(x)
  if (is.null(weights))   weights <- rep.int(1, nobs)
  if (is.null(offset))    offset <- rep.int(0, nobs)
  newfit = .Call("glm_bas",
    			RX=x, 
			RY = y,
    			family=family, 
			Roffset = offset,
    			Rweights = weights, 
			Rcontrol=control, 
			PACKAGE="BAS")
  newfit$fitted.values <- newfit$mu
  return(newfit)
}
#his this function for now
.gglm.lpy <- function(X,y,a,b,s,fit = NULL) {
	n = dim(X)[1];
  	p = dim(X)[2];
	if (is.null(fit)) { # for now
		fit <- gglm.fit(cbind(1,X),y);
	}
	if (is.null(dim(X))) {
		dim(X) <- c(length(y),1)
	}
	lpY <-  .Call("gglm_lpy",
    			RX=X, 
			RY = y,
    			Ra = a, 
			Rb = b,
    			Rs = s, 
			Rcoef = fit$coef, 
			Rmu = fit$mu,
			PACKAGE="BAS")
 }

#to be reomved later
.lpy.gglm_old <- function(X,y,a,b,s,fit = NULL) {
	n = dim(X)[1];
  	p = dim(X)[2];
	if (is.null(fit)) { # for now
		if (p == 0) {
			fit = glm(y ~ 1, family = binomial(link = 'logit'));
		} else {
			fit = glm(y ~ X, family = binomial(link = 'logit'));
		}
	}
 	lpY <- .lpY.logit(X, y, a, b, s, fit);
}

