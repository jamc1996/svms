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
}

void setH(struct Fullproblem *prob, struct denseData *ds)
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
    std::cout << "Alpha hat = " << sp->alphaHat[i] << '\n';
  }

  for (int i = 0; i < sp->p; i++) {
    sp->rho[i] = sp->alphaHat[i];
    for (int j = 0; j < sp->p; j++) {
      sp->rho[i] -= ((sp->yHat[i]*sp->yHat[j])/((double)(sp->p)))*sp->alphaHat[i];
    }
  }

  for (int i = 0; i < sp->p; i++) {
    std::cout << "pal = " << sp->rho[i] << '\n';
  }

  for (int i = 0; i < sp->p; i++) {
    fp->alpha[fp->active[i]] = sp->rho[i];
  }

  for (int i = 0; i < fp-> q; i++) {
    for (int j = 0; j < sp->p; j++) {
      fp->gradF[fp->inactive[i]] -= fp->partialH[i][j] * sp->rho[j];
    }
  }

  for (int i = 0; i < fp->p; i++) {
    for (int j = 0; j < fp->p; j++) {
      fp->gradF[fp->active[i]] = sp->gamma[i];
    }
  }

}

int partialSwap(struct Fullproblem *fp, struct Projected *sp)
{
  int badVals = 0;
  for (int i = 0; i < sp->p; i++) {
    if (fabs(sp->alphaHat[i]-sp->C) < 0.00001 || sp->alphaHat[i] < 0.000001) {
      badVals++;
    }
  }
  int* pLocations = (int*)malloc(sizeof(int)*badVals);
  int* qLocations = (int*)malloc(sizeof(int)*badVals);
  double* betaVal = (double*)malloc(sizeof(double)*badVals);
  int* qIndex = (int*)malloc(sizeof(int)*badVals);
  int j = 0;
  for (int i = 0; i < sp->p; i++) {
    if (fabs(sp->alphaHat[i]-sp->C) < 0.00001 || sp->alphaHat[i] < 0.000001) {
      pLocations[j] = i;
      qLocations[i] = -1;
      betaVal[i] = 100.0;
    }
  }
  for (int i = 0; i < fp->q; i++) {
    for (int j = 0; j < badVals; j++) {
      if (fp->beta[i]<betaVal[j]) {
        for (int l = badVals - 1; l > j; l--) {
          qLocations[l] = qLocations[l-1];
          betaVal[l] = betaVal[l-1];
        }
        qLocations[j] = i;
        betaVal[j] = fp->beta[i];
        break;
      }
    }
  }
  for (int i = 0; i < badVals; i++) {
std::cout << "/* message */ " << qLocations[i] << '\n';
}

  for (int i = 0; i < badVals; i++) {
    qIndex[i] = fp->inactive[qLocations[i]];
    fp->inactive[qLocations[i]] = fp->active[pLocations[i]];
    fp->active[pLocations[i]] = qIndex[i];
  }

  free(qLocations);
  free(pLocations);
  free(qIndex);
  free(betaVal);

  return 1;

}

int swapMostNegative(struct Fullproblem *fp)
{
  int* location = (int*)malloc(sizeof(int)*fp->p);
  int* index = (int*)malloc(sizeof(int)*fp->p);
  double* betaVal = (double*)malloc(sizeof(double)*fp->p);
  for (int i = 0; i < fp->p; i++) {
    location[i] = -1;
    betaVal[i] = 100.0;
  }
  for (int i = 0; i < fp->q; i++) {
    for (int j = 0; j < fp->p; j++) {
      if (fp->beta[i]<betaVal[j]) {
        for (int k = (fp->p) - 1; k > j ; k--) {
          location[k] = location[k-1];
          betaVal[k] = betaVal[k-1];
        }
        location[j] = i;
        betaVal[j] = fp->beta[i];
        break;
      }
    }
  }
  if (betaVal[0]>=0.0) {
    free(location);
    free(index);
    free(betaVal);
    return 0;
  }

  for (int i = 0; i < fp->p; i++) {
    index[i] = fp->inactive[location[i]];
    fp->inactive[location[i]] = fp->active[i];
    fp->active[i] = index[i];
  }


  free(location);
  free(index);
  free(betaVal);

  return 1;
}


void calculateBeta(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds)
{
  for (int i = 0; i < fp->q; i++) {
    fp->beta[i] = fp->gradF[fp->inactive[i]] - sp->ytr*ds->y[fp->inactive[i]];
  }
}
