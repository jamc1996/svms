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

  // Input processed:
  parse_arguments(argc, argv, &filename, &parameters);
  read_file(filename, &ds);
  preprocess(&ds);

  // Projected problem size chosen temporarily
  int p = 8;

  // Full problem allocated and filled in:
  alloc_prob(&fullP, ds.nInstances, p);
  init_prob(&fullP, ds.nInstances, p, &ds);

  // Subproblem allocated:
  alloc_subprob(&subP, p, &fullP, &ds);

  // We loop until no negative entries in beta:
  int k = 1;
  int max_iters = 10;
  int i = 0;
  while(k){
    init_subprob(&subP, &fullP, &ds);
    cg(&subP);
    calcYTR(&subP);
    updateAlphaR(&fullP, &subP);
    calculateBeta(&fullP, &subP, &ds);
    std::cout << "chaning::::" << '\n';
    for (int i = 0; i < fullP.p; i++) {
      std::cout << fullP.active[i] << '\n';
    }

    swapMostNegative(&fullP);
    updatePartialH(&fullP, &ds);
    
    for (int i = 0; i < fullP.n; i++) {
      std::cout << "alpha =" << fullP.alpha[i] << '\n';
    }
    std::cout << '\n';

    i++;
    if(i == max_iters){
      break;
    }

  }

  //Memory freed
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
