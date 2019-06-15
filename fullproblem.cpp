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

  for (int i = 0; i < p/2; i++) {
    prob->active[i] = i;
  }

  for (int i = 0; i < p/2; i++) {
    prob->active[i+p/2] = ds->nPos + i;
  }

  if (p%2 != 0) {
    prob->active[p-1] = ds->nPos+p/2;
    std::cout << "fullproblem.cpp: init_prob(): Not yet working" << '\n';
    exit(1);
  }
  else{
    for (int i = p/2; i < ds->nPos; i++) {
      prob->inactive[i-(p/2)] = i;
    }
    for (int i = p/2; i < ds->nNeg; i++) {
      prob->inactive[ds->nPos-(p)+i] = ds->nPos + i;
    }
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
    std::cout << "fp acti s " << fp->active[i] << '\n';
  }

  for (int i = 0; i < sp->p; i++) {
    fp->alpha[fp->active[i]] = sp->rho[i];
  }

  for (int i = 0; i < fp-> q; i++) {
    std::cout << fp->inactive[i] << '\n';
    for (int j = 0; j < sp->p; j++) {
      fp->gradF[fp->inactive[i]] -= fp->partialH[i][j] * sp->rho[j];
    }
  }

  for (int i = 0; i < fp->p; i++) {

    for (int j = 0; j < i; j++) {
      fp->gradF[fp->active[i]] -= sp->H[j][i]*sp->rho[j];
    }
    for (int j = i; j < fp->p; j++) {
      fp->gradF[fp->active[i]] -= sp->H[i][j]*sp->rho[j];
    }
  }

}

int singleswap(struct denseData *ds, struct Fullproblem *fp, struct Projected *sp, int n)
{
  int worst = -1;
  int target, change=1;
  double tester = 100.0;

  if (n >= sp->p) {
    change = -1;
    n-=sp->p;
  }
  target = ds->y[fp->active[n]]*change;
  std::cout << "n is " << n << '\n';
  std::cout << "change is " << change << '\n';
  std::cout << "target is " << target << '\n';

  target = change*ds->y[fp->active[n]];

  for (int i = 0; i < fp->q; i++) {
    if( ds->y[fp->inactive[i]] == target )  {
      if (fp->beta[i]<tester) {
        worst = i;
        tester = fp->beta[i];
      }
    }
  }

  std::cout << "worst = " << worst << '\n';
  int temp = fp->active[n];


  if (change < 0)
  {
    adjustGradF(fp, ds, n, worst, change);
    fp->alpha[fp->inactive[worst]]+= sp->C - fp->alpha[fp->active[n]];
    fp->alpha[fp->active[n]] = sp->C ;
    fp->active[n] = fp->inactive[worst];
    fp->inactive[worst] = temp;
  }
  else
  {
    fp->alpha[fp->inactive[worst]] += fp->alpha[fp->active[n]];
    fp->alpha[fp->active[n]] = 0;
    fp->active[n] = fp->inactive[worst];
    fp->inactive[worst] = temp;
  }


  return 1;

}

void adjustGradF(struct Fullproblem *fp, struct denseData *ds, int n, int worst, int signal)
{
  if (signal == -1) {
    for (int i = 0; i < fp->q; i++) {
      fp->gradF[fp->inactive[i]] += fp->partialH[i][n]*fp->alpha[n] ;
    }
  }
  else
  {
    for (int i = 0; i < fp->q; i++) {
      fp->gradF[fp->inactive[i]] -= fp->partialH[i][n]*fp->alpha[n] ;
    }
  }
  double* temp = (double*) malloc(sizeof(double)*fp->q);
  for (int j = 0; j < fp->q; j++) {
    temp[j] = 0.0;
    for (int k = 0; k < ds->nFeatures; k++) {
      temp[j]+=ds->data[fp->inactive[worst]][k]*ds->data[fp->inactive[j]][k];
    }
    temp[j]*=ds->y[fp->inactive[worst]]*ds->y[fp->inactive[j]];
  }
  free(temp);
}

int swapMostNegative(struct Fullproblem *fp)
/*  If the cg has fully converged we swap all p values out for a new set
 *  based on how negative their beta value is.    */
{
  int* location = (int*)malloc(sizeof(int)*fp->p);
  int* index = (int*)malloc(sizeof(int)*fp->p);
  double* betaVal = (double*)malloc(sizeof(double)*fp->p);
  for (int i = 0; i < fp->p; i++)
  {
    location[i] = -1;
    betaVal[i] = 100.0;
  }
  for (int i = 0; i < fp->q; i++)
  {
    for (int j = 0; j < fp->p; j++)
    {
      if (fp->beta[i]<betaVal[j])
      {
        for (int k = (fp->p) - 1; k > j ; k--)
        {
          location[k] = location[k-1];
          betaVal[k] = betaVal[k-1];
        }
        location[j] = i;
        betaVal[j] = fp->beta[i];
        break;
      }
    }
  }
  if (betaVal[0]>=0.0)
  {
    free(location);
    free(index);
    free(betaVal);
    return 0;
  }

  for (int i = 0; i < fp->p; i++)
  {
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
/* Function to calculate the beta vector - tracks how suboptimal the values for*/
{
  for (int i = 0; i < fp->q; i++) {
    fp->beta[i] = fp->gradF[fp->inactive[i]] - (sp->ytr*ds->y[fp->inactive[i]]/(double)(fp->p));
  }
}
