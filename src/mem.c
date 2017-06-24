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

long *lvecalloc(int nr)
{
  long *x;
  x=(long *) R_alloc((unsigned) nr,sizeof(long));
  return x;
}



float *fvecalloc(int nr)
{
  float *x;
  x=(float *) R_alloc((unsigned) nr,sizeof(float));
  return x;
}

int *ivecalloc(int nr)
{
  int *x;
  x= (int *) R_alloc( nr, sizeof(int));
  return x;
}

/*
double **matallocold(int nr, int nc)
{
int k;
double **x;
x= Calloc( nr, double *);
for (k=0;k<nr;k++){
  x[k]= Calloc(nc, double);
  }
return x;
}
*/

float **fmatalloc(int nr, int nc)
{
int k;
float **x;
x=(float **) R_alloc((unsigned) nr, sizeof(float *));

for (k=0;k<nr;k++){
  x[k]=(float *) R_alloc((unsigned) nc, sizeof(float));
}
return x;
}

double **matalloc(int nr, int nc)
{
int k;
double **x;
x=(double **) R_alloc((unsigned) nr, sizeof(double *));

for (k=0;k<nr;k++){
  x[k]=(double *) R_alloc((unsigned) nc, sizeof(double));
}
return x;
}

int  **imatalloc(int nr, int nc)
{
int k;
int **x;
x = (int **) R_alloc(nr, sizeof(int *));
for (k=0;k<nr;k++){
  x[k] = (int *) R_alloc(nc, sizeof(int));
}
return x;
}

unsigned char  **cmatalloc(int nr, int nc)
{
int k;
unsigned char  **x;
x = (unsigned char **) R_alloc(nr, sizeof(unsigned char *));
for (k=0;k<nr;k++){
  x[k] = (unsigned char  *) R_alloc(nc, sizeof(unsigned char));
}
return x;
}


void  freemat(double **mat, int  nr)
{
int k;

for (k=0;k<nr;k++){
  Free(mat[k]);
  }
Free(mat);

}

void  freeimat(int **mat, int  nr)
{
int k;

for (k=0;k<nr;k++){
  Free(mat[k]);
  }
Free(mat);

}

void  freechmat(unsigned char **mat, int  nr)
{
int k;

for (k=0;k<nr;k++){
  Free(mat[k]);
  }
Free(mat);
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

