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
    prob->inactive[i-p] = i;
  }

  for (int i = 0; i < prob->q; i++) {
    for (int j = 0; j < prob->p; j++) {
      prob->partialH[i][j] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        prob->partialH[i][j]+=ds->data[prob->inactive[i]][k]*ds->data[prob->active[j]][k];
      }
      prob->partialH[i][j]*=ds->y[prob->inactive[i]]*ds->y[prob->active[j]];
    }
  }
}

void updatePartialH(struct Fullproblem *prob, struct denseData *ds)
{
  for (int i = 0; i < prob->q; i++) {
    for (int j = 0; j < prob->p; j++) {
      prob->partialH[i][j] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        prob->partialH[i][j]+=ds->data[prob->inactive[i]][k]*ds->data[prob->active[j]][k];
      }
      prob->partialH[i][j]*=ds->y[prob->inactive[i]]*ds->y[prob->active[j]];
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

  for (int i = 0; i < fp-> q; i++) {
    for (int j = 0; j < sp->p; j++) {
      fp->gradF[fp->inactive[i]] -= fp->partialH[i][j] * sp->rho[j];
    }
  }

}

int swapMostNegative(struct Fullproblem *fp)
{
  int* temp = (int*)malloc(sizeof(int)*fp->p);
  int* temp2 = (int*)malloc(sizeof(int)*fp->p);
  double* tempD = (double*)malloc(sizeof(double)*fp->p);
  for (int i = 0; i < fp->p; i++) {
    temp[i] = -1;
    tempD[i] = 100.0;
  }
  for (int i = 0; i < fp->q; i++) {
    for (int j = 0; j < fp->p; j++) {
      if (fp->beta[i]<tempD[j]) {
        for (int k = (fp->p) - 1; k > j ; k--) {
          temp[k] = temp[k-1];
          tempD[k] = tempD[k-1];
        }
        temp[j] = i;
        tempD[j] = fp->beta[i];
        break;
      }
    }
  }

  for (int i = 0; i < fp->p; i++) {
    temp2[i] = fp->inactive[temp[i]];
    fp->inactive[temp[i]] = fp->active[i];
    fp->active[i] = temp2[i];
  }

  free(temp);
  free(temp2);
  free(tempD);
  if (fp->active[0]>=0) {
    return 0;
  }
  return 1;
}


void calculateBeta(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds)
{
  for (int i = 0; i < fp->q; i++) {
    fp->beta[i] = fp->gradF[fp->inactive[i]] - sp->ytr*ds->y[fp->inactive[i]];
  }
}
