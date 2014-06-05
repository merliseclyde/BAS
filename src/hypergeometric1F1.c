extern double hyperg(double, double, double);

void hypergeometric1F1(double *a, double *b, double *x, double *y, int *npara)
{ 
  int k;
  for (k = 0; k < *npara; k++) {
    y[k] = hyperg(a[k], b[k], x[k]);
  }
}

 
