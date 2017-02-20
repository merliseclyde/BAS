extern double hyp2f1(double, double, double, double);

// [[register]]
void hypergeometric2F1(a, b, c, x, y)
     double *y, *x, *a, *b, *c;
{

 *y = hyp2f1(*a, *b, *c, *x);
}

 
