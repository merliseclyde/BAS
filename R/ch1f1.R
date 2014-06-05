.ch1f1 = function(a, b, x) {
  #rpois(1, 1);
  if(length(a)!= length(b) || length(a)!= length(x)){
  	cat('Error in "ch1f1": lengths of vector a, b, x are not equal!');
  }

  npara = length(a);
  
  ## Since Linex system tends to report error for negative x, we use the following fomular to convert it to positive value
  ## 1F1(a, b, x) = 1F1(b - a, b, -x) * exp(x) 
  negativex = which(x < 0);
  a[negativex] = b[negativex] - a[negativex];
  x[negativex] = -x[negativex];
  
  results.1f1 = .C('hypergeometric1F1', as.double(a), as.double(b), as.double(x), double(npara), as.integer(npara), PACKAGE = 'BAS')[[4]];
  results.1f1[negativex] = results.1f1[negativex] * exp(-x[negativex]);
  
  return( results.1f1 ); 
}
