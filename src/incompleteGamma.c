extern double igam(double, double);

void incompleteGamma(double *a, double *x, double *y)
{
 *y = igam(*a, *x);
}

 
