/* Modified to be compatible with R memory allocation */
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "bas.h"

double *vecalloc(int nr)
{
  double *x;
  x=(double *) R_alloc(nr,sizeof(double));
  return x;
}


int *ivecalloc(int nr)
{
  int *x;
  x= (int *) R_alloc( nr, sizeof(int));
  return x;
}




unsigned char  **cmatalloc(int nr, int nc)
{
int k;
unsigned char  **x;
x = (unsigned char **) R_alloc(nr, sizeof(unsigned char *));
for (k=0;k<nr;k++){
  x[k] = (unsigned char  *) R_alloc(nc, sizeof(unsigned char));
  memset(x[k], 0, nc*sizeof(unsigned char));
}
return x;
}



/* getListElement from Writing R Extensions */


     SEXP getListElement(SEXP list, char *str)
     {
       SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
       int i;

       for (i = 0; i < length(list); i++)
         if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
         }
       return elmt;
     }


