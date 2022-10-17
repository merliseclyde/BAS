extern double hyp2f1(double, double, double, double);


void hypergeometric2F1(double *a, double  *b, double *c, double *x, double *y)
     {

 *y = hyp2f1(*a, *b, *c, *x);
}

 
