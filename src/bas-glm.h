void insert_model_tree(struct Node *tree, struct Var *vars,  int n, int *model, int num_models);

int *GetModel_m(SEXP Rmodel_m, int *model, int p); 
	
void CreateTree(NODEPTR branch, struct Var *vars, int *bestmodel, int *model, int n, int m, SEXP modeldim);

void CreateTree_with_pigamma(NODEPTR branch, struct Var *vars, int *bestmodel, int *model, int n, int m, 
							 SEXP modeldim, double *pigamma);

double GetNextModelCandidate(int pmodel_old, int n, int n_sure, int *model, struct Var *vars, double problocal,
							 int *varin, int *varout);

void GetNextModel_swop(NODEPTR branch, struct Var *vars, int *model, int n, int m,  double *pigamma, 
					   double problocal, SEXP modeldim,int *bestmodel);

void Substract_visited_probability_mass(NODEPTR branch, struct Var *vars, int *model, int n, int m, double *pigamma, double eps);

void SetModel1(SEXP Rfit, SEXP Rmodel_m, 
	       SEXP beta, SEXP se, SEXP modelspace, SEXP deviance, SEXP R2, SEXP Q, SEXP Rintercept, int m);

void SetModel2(double logmargy, double shrinkage_m, double prior_m,
			  SEXP sampleprobs, SEXP logmarg, SEXP shrinkage, SEXP priorprobs, int m);

SEXP glm_FitModel(SEXP RX, SEXP RY, SEXP Rmodel_m,  //input data
		  SEXP Roffset, SEXP Rweights,
		  glmstptr * glmfamily, SEXP Rcontrol,
		  SEXP Rlaplace, betapriorptr * betapriorfamily);

SEXP glm_bas(SEXP RX, SEXP RY, glmstptr * family, SEXP Roffset, SEXP Rweights, SEXP Rcontrol);

SEXP gglm_lpy(SEXP RX, SEXP RY,SEXP Rcoef, SEXP Rmu, glmstptr * glmfamily, betapriorptr * betapriorfamily, SEXP Rlaplace);

