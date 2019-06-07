#include "fullproblem.h"
#include <iostream>

void alloc_prob(struct Fullproblem *prob, int n, int p)
{
  prob->n = n;
  prob->p = p;
  prob->q = n-p;

  prob->alpha = (double*)calloc(prob->n, sizeof(double) );
  prob->gradF = (double*)malloc(sizeof(double)*n);

  prob->active = (int*)malloc(sizeof(int)*p);
  prob->inactive = (int*)malloc(sizeof(int)*prob->q);



  prob->beta = (double*)malloc(sizeof(double)*(prob->q));



  prob->h = (double*)malloc(sizeof(double)*prob->q*prob->p);
  prob->partialH = (double**)malloc(sizeof(double*)*prob->q);
  for (int i = 0; i < prob->q; i++) {
    prob->partialH[i] = &prob->h[i*prob->p];
  }

}

void init_prob(struct Fullproblem *prob, int n, int p, denseData *ds)
{
  for (int i = 0; i < n; i++) {
    prob->gradF[i] = 1.0;
  }

  for (int i = 0; i < p; i++) {
    prob->active[i] = i;
  }
  for (int i = p; i < n; i++) {
    prob->inactive[i] = i;
  }

  for (int i = 0; i < prob->q; i++) {
    for (int j = 0; j < prob->p; j++) {
      prob->partialH[i][j] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        prob->partialH[i][j]+=ds->data[i+p][k]*ds->data[j][k];
      }
    }
  }
}

void updateAlphaR(struct Fullproblem *fp, struct Projected *sp)
{
  
  for (int i = 0; i < sp->p; i++) {
    sp->rho[i] = sp->alphaHat[i];
    for (int j = 0; j < sp->p; j++) {
      sp->rho[i] -= (sp->yHat[i]*sp->yHat[j])/((double)(sp->p))*sp->alphaHat[i];
    }
  }

  for (int i = 0; i < sp->p; i++) {
    fp->alpha[fp->active[i]] += sp->rho[i];
  }

  for (int i = 0; i < fp->n; i++) {
    for (int j = 0; j < sp->p; j++) {
      fp->gradF[i] -= fp->partialH[i][j] * sp->rho[j];
    }
  }
}

void calculateBeta(struct Fullproblem *fp, struct Projected *sp)
{
  for (int i = 0; i < fp->q; i++) {
    fp->beta[i] -= sp->ytr;
  }
}
