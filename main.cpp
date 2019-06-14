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
  int max_iters = 2;
  int itt = 0;
  int n = 0;

  while(k){

    // H matrix columns fixed and subproblem changed
    setH(&fullP, &ds);
    init_subprob(&subP, &fullP, &ds);

    // Conjugate gradient algorithm with interrupt if
    n = cg(&subP);

    // Computations for choosing next elements to analyze
    calcYTR(&subP);
    updateAlphaR(&fullP, &subP);
    calculateBeta(&fullP, &subP, &ds);
std::cout << "n is " << n << '\n';
    if (n) {
      // Case where interrupt needed:
      k = singleswap(&ds, &fullP,&subP,n);
      for (int i = 0; i < subP.p; i++) {
        std::cout << "activia " << fullP.active[i] << '\n';
      }
    }
    else{
      //  Without interrupt we swap in full
      k = swapMostNegative(&fullP);
      for (int i = 0; i < subP.p; i++) {
        std::cout << "muller " << fullP.active[i] << '\n';
      }
    }

    for (int k = 0; k < fullP.n; k++) {
      std::cout << " alpha [ " << k << " ] =  " << fullP.alpha[k] << '\n';
    }

    itt++;
    if(itt == max_iters){
      break;
    }
  }
  if (k==0) {
    std::cout << "k goes to 0." << '\n';
  }
  //Memory freed
  freeDenseData(&ds);
  freeFullProblem(&fullP);
  freeSubProblem(&subP);

  return 0;
}

void   freeSubProblem( struct Projected* subP)
/* Function to free dynamically allocated memory in subproblem stuct. */
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
/* Function to free dynamically allocated memory in dense data set struct. */
{
  free(ds->data);
  free(ds->data1d);
  free(ds->y);
  free(ds->instanceLabels);
  free(ds->featureLabels);
}

void  freeFullProblem(struct Fullproblem *fp)
/* Function to free dynamically allocated memory in fullproblem struct */
{
  free(fp->alpha);
  free(fp->beta);
  free(fp->gradF);

  free(fp->active);
  free(fp->inactive);

  free(fp->partialH);
  free(fp->h);
}
