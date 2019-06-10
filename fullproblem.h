#ifndef FULLPROBLEM_H
#define FULLPROBLEM_H

#include <stdlib.h>

#include "svm.h"

void alloc_prob(struct Fullproblem *prob, int n, int p);
void init_prob(struct Fullproblem *prob, int n, int p, denseData *ds);
void updateAlphaR(struct Fullproblem *fp, struct Projected *sp);
void calculateBeta(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds);
int swapMostNegative(struct Fullproblem *fp);
void updatePartialH(struct Fullproblem *prob, struct denseData *ds);

#endif
