#########################################################################
###### For CHg prior: for any (sub)model M
###### 1) Calculate approximate marginal likelihood
###### 2) Here p and X are sub-model's, not full model's!!!
#########################################################################

###### Inputs ######
## a, b, s: scale
###### Outputs ######
## lpY: marginal likelihood, for a vector with the same length to (a, b, s); contains infomration from a, b, s
.lpY.logit = function(X, y, a, b, s, fit){		
  	n = dim(X)[1];
  	p = dim(X)[2];
  	if(length(a)!= 1 || length(b)!= 1 || length(s)!= 1){
  		cat('Error in "lpY.logit": a,b,s must be scalers');
  	}
  	lpY = NA;     
  
  	## if null model
  	if(p == 0){
		mu <- fit$fitted.values; # = odds ratio = mu
	
    		## log maximized likelihood
		lC <- sum( dbinom(y, 1, mu, log = TRUE));
		Ieta <-  mu * (1 - mu);
		sum.Ieta <- sum(Ieta);
		lpY = lC + 0.5 * log(2 * pi) - 0.5 * log(sum.Ieta);    
    		## for Jeffreys's prior (b = 0), set lpY to NA, since then CH g-prior doesn't appropriate for the null model in this case
    		if(b == 0){
      		lpY = NA;
    		}
  	}
  
  	## if not null model	
  	if(p > 0){
    		## mle estimators
    		beta = fit$coef[-1];
		mu <- fit$fitted.values; # = odds ratio = mu
		exp_eta <- mu / (1-mu);

    		## log maximized likelihood
		lC <- sum( dbinom(y, 1, mu, log = TRUE));

    		## "centering"
		Ieta <-  mu * (1 - mu);
		sum.Ieta <- sum(Ieta);
		Xc <- X - rep(1,n) %*% t((t(X) %*% Ieta)) / sum.Ieta;
		Q <- sum((Xc %*% beta)^2 * Ieta);		
		#invIbeta <- solve(t(Xc) %*% sweep(Xc,1,Ieta,FUN = "*"));

    		## marginal likelihood
    		lpY = lC + 0.5 * log(2 * pi) - 0.5 * log(sum.Ieta) + lbeta((b + p) / 2, a / 2) + log(.ch1f1((b + p) / 2, (a + b + p) / 2, -  (s + Q) / 2))  - log(.ch1f1(b / 2, (a + b) / 2, - s / 2));
    		## doesn't apply for the Jeffreys prior
		if (a > 0 & b > 0) {
			lpY = lpY - lbeta(b / 2, a / 2);
		}
    		
  	}
  	return(lpY);
}
