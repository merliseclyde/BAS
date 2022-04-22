v = 1; s = 0; k = 0; 


gordy = function(u, a = 3, b = 1, r = 0, s = 1, v = 1, k = 1) {
   v^(a/2)* exp(s/v)* u^{a/2 - 1}*(1 - v*u)^{b/2 - 1}*(k + (1 - k)*u*v)^{-r}*exp(-s*u)
}

integrate(gordy, 0 , 1)
phi1_dens = function(u, a,b,c,x,y) {u^{a-1}*(1 - u)^{c-a-1}*(1 - x*y)^{-b}*exp(y*u)}

phi1_int = function(a,b,c,x=1,y=0) {
   (gamma(c)/gamma(a)*gamma(c-a)) *
    integrate(function(u){u^{a-1}*(1 - u)^{c-a-1}*(1 - y*u)^{-b}*exp(x*u)}, 0, 1)
}

hyper.gn.dens = function(g, p, n, Q, a = 4, s = 0) {
  (1 + g/n)^{-a/2} * ( 1 + g)^{-p/2} * exp(-(s + Q)/2)
}

robust.dens = function(g, p, n, Q) {
  v= (n+1)/(p+1)
  u = 1/(1 + g)
  s = 0; a = 3; b = 0; r = 1.5; k = 1
  (1+g)^{-2}*u^{(a + p)/2 - 1}*(1 - v*u)^{b/2 - 1}*(k + (1 - k)*u*v)^{-r}*exp(-(s + Q)*u)*(u < 1/v)
}

ttch.dens = function(g, p, n, Q, r = 1, b = 1, k = 1, v=NULL) {
   if (is.null(v))  v= (n+1)/(p+1)
  u = 1/(1 + g)
  s = 0; a = 3;
 u^{a/2 + p/2 - 1}*(1 - v*u)^{b/2 - 1}*(k + (1 - k)*u*v)^{-r}*exp(-(s + Q)*u)*(u < 1/v)
}

p = 97
153.8388
Q = 55.51288
-114.5882