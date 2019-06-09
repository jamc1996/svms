#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"

#include <iostream>

void freeDenseData(struct denseData *ds);
void  freeFullProblem(struct Fullproblem *fp);
void   freeSubProblem( struct Projected* subP);

int main(int argc, char *argv[]) {
  char* filename = NULL;
  struct svm_args parameters;
  struct denseData ds;
  struct Fullproblem fullP;
  struct Projected subP;

  parse_arguments(argc, argv, &filename, &parameters);
  read_file(filename, &ds);
  preprocess(&ds);

  for (int i = 0; i < ds.nInstances; i++) {
    for (int j = 0; j < ds.nFeatures; j++) {
    //  std::cout << ds.data[i][j] << '\t';
    }
    //std::cout << ds.y[i] << '\n';
  }

  int p = 8;

  alloc_prob(&fullP, ds.nInstances, p);
  init_prob(&fullP, ds.nInstances, p, &ds);

  for (int i = 0; i < fullP.q; i++) {
    for (int j = 0; j < fullP.p; j++) {
      std::cout << fullP.partialH[i][j] << '\t';
    }
    std::cout << '\n';
  }

  alloc_subprob(&subP, p, &fullP, &ds);
  init_subprob(&subP, &fullP, &ds);

  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
      std::cout << subP.H[i][j] << '\t';
    }
    std::cout << '\n';
  }

  cg(&subP);
  for (int i = 0; i < p; i++) {
    std::cout << subP.alphaHat[i] << '\n';
  }
  calcYTR(&subP);


  updateAlphaR(&fullP, &subP);
  calculateBeta(&fullP, &subP, &ds);


  for (int i = 0; i < fullP.q; i++) {
    std::cout << "beta is " << fullP.beta[i] << '\n';
  }
  for (int i = 0; i < fullP.n; i++) {
    //std::cout << "gradF is " << fullP.gradF[i] << '\n';
  }

  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
  //    std::cout << subP.H[i][j] << '\t';
    }
    std::cout << '\n';
  }

  swapMostNegative(&fullP);


  freeDenseData(&ds);
  freeFullProblem(&fullP);
  freeSubProblem(&subP);

  return 0;
}

void   freeSubProblem( struct Projected* subP)
{
  free(subP->alphaHat);
  free(subP->yHat);
  free(subP->rHat);
  free(subP->H);

  free(subP->gamma);
  free(subP->rho);
  free(subP->Hrho);

  free(subP->h);
}

void freeDenseData(struct denseData *ds)
{
  free(ds->data);
  free(ds->data1d);
  free(ds->y);
  free(ds->instanceLabels);
  free(ds->featureLabels);
}

void  freeFullProblem(struct Fullproblem *fp)
{
  free(fp->alpha);
  free(fp->beta);
  free(fp->gradF);

  free(fp->active);
  free(fp->inactive);

  free(fp->partialH);
  free(fp->h);
}
