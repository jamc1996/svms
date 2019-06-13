#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include "svm.h"

#include <stdlib.h>
#include <iostream>

void alloc_subprob(struct Projected *sp, int p, struct Fullproblem *fp, struct denseData *ds);
void init_subprob(struct Projected *sp, struct Fullproblem *fp, struct denseData *ds);
void init_symmetric(struct Projected *sp, int p);

int cg(struct Projected *sp);

void init_error(struct Projected* sp);
void calc_Hrho(struct Projected *sp);
void linearOp2(double* rho, double* gama, double mu, int p);

void linearOp(double* alpha, double* rho, double lambda, int p);
void updateGamma(struct Projected *sp, double lambda);
void calcYTR(struct Projected *sp);
void copy_vector(double* a, double* b, int p);

void constraint_projection(double* vecOut, double* vecIn, double* y, int p);

double inner_prod(double *a, double *b, int p);

#endif
