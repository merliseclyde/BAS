#include "bas.h"

double GetNextModelCandidate(int pmodel_old, int n, int n_sure, int *model, struct Var *vars, double problocal,
                             int *varin, int *varout, SEXP Rparents) {
  double MH = 1.0;
  if (pmodel_old == n_sure || pmodel_old == n_sure + n){
    MH =  random_walk_heredity(model, vars,  n, Rparents);
    MH =  1.0 - problocal;
  } else {
    if (unif_rand() < problocal) {
      // random
      MH =  random_switch_heredity(model, vars, n, pmodel_old, varin, varout, Rparents );
    } else {
      // Random walk proposal flip bit//
      MH =  random_walk_heredity(model, vars,  n, Rparents);
    }
  }
  return MH;
}


// # nocov start
// code below not being used internally - legacy?

double random_walk(int *model, struct Var *vars, int n) {
  int index;
  index = ftrunc(n*unif_rand());
  model[vars[index].index] = 1 - model[vars[index].index];
  return(1.0);
}

double random_switch(int *model, struct Var *vars, int n, int pmodel, int *varin, int *varout) {
  int  j, k, swapin, swapout, num_to_swap_in, num_to_swap_out;
  
  
  j = 0; k = 0;
  while (j < n && k < pmodel)
  {
    if (model[vars[j].index]==1) {varin[k] = vars[j].index; k++;}
    j++ ;
  }
  num_to_swap_in = k;
  j = 0; k = 0;
  
  while (j< n)
  {
    if (model[vars[j].index]==0) {varout[k] = vars[j].index; k++;}
    j++ ;
  }
  num_to_swap_out = k;
  
  swapin = ftrunc(unif_rand()*num_to_swap_in);    // swapin :corresponds to position of randomly chosen included variable
  swapout = ftrunc(unif_rand()*num_to_swap_out);  // swapout :corresponds to position of randomly chosen excluded variable
  
  model[varin[swapin]] = 0;
  model[varout[swapout]] =1;
  
  
  return(1.0);
}


// # nocov end

double random_walk_heredity(int *model, struct Var *vars, int n, SEXP Rparents) {
  int index,p,j;
  double *parents;
  
  parents = REAL(Rparents);
  
  index = ftrunc(n*unif_rand());
  model[vars[index].index] = 1 - model[vars[index].index];
  
  
  int *dims = INTEGER(getAttrib(Rparents,R_DimSymbol));
  p = dims[0];
  
  
  if (p > 1) {
    // force in parents
    //  Rprintf("%d %d %d %d\n",n,p,  vars[index].index, model[vars[index].index]);
    if (model[vars[index].index] == 1) {
      //  traverse row index of parents to add any missing parents/sibs
      for (j = 0; j < p; j++) {
        //     Rprintf("%d ", (int) parents[vars[index].index*p + j]);
        if (parents[vars[index].index + p*j] == 1.0) {
          model[j] =  model[vars[index].index];
        }
      }}
    else {
      //  to drop index, traverse column of parents to identify children/sibs
      //  that also need to be dropped.
      for (j = 0; j < p; j++) {
        //      Rprintf("%d ", (int) parents[vars[index].index + p*j]);
        if (parents[vars[index].index*p +j] == 1.0) {
          model[j] =  model[vars[index].index];
        } }
    }
    //    Rprintf("\n");
    
  }
  return(1.0);
}

double random_switch_heredity(int *model, struct Var *vars, int n,
                              int pmodel, int *varin, int *varout, SEXP Rparents)
{
  int  j, k, p, swapin, swapout, num_to_swap_in, num_to_swap_out;
  double *parents;
  
  j = 0; k = 0;
  while (j < n && k < pmodel)
  {
    if (model[vars[j].index]==1) {
      varin[k] = vars[j].index;
      k++;}
    j++ ;
  }
  num_to_swap_in = k;
  j = 0; k = 0;
  
  while (j< n)
  {
    if (model[vars[j].index]==0) {
      varout[k] = vars[j].index;
      k++;}
    j++ ;
  }
  num_to_swap_out = k;
  
  swapin = ftrunc(unif_rand()*num_to_swap_in);    // swapin :corresponds to position of randomly chosen included variable
  swapout = ftrunc(unif_rand()*num_to_swap_out);  // swapout :corresponds to position of randomly chosen excluded variable
  
  model[varin[swapin]] = 0;
  model[varout[swapout]] = 1;
  
  parents = REAL(Rparents);
  int *dims = INTEGER(getAttrib(Rparents,R_DimSymbol));
  p = dims[0];
  
  // force in parents and sibs of variable that was swapped in
  
  if (p > 1) {
    //  to drop swapin, traverse column of parents to identify children/sibs
    //  that also need to be dropped.  ignore others
    for (j = 0; j < p; j++) {
      if (parents[varin[swapin]*p +j] == 1.0)   model[j] = 0;
    }
    
    //  now traverse row of added variable in parents to add any missing parents/sibs
    for (j = 0; j < p; j++) {
      if (parents[varout[swapout] + p*j] == 1.0)   model[j] = 1;
    }
  }
  return(1.0);
}

