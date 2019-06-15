#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"

#include <iostream>

void freeDenseData(struct denseData *ds);
void  freeFullproblem(struct Fullproblem *fp);
void   freeSubProblem( struct Projected* sp);

int main(int argc, char *argv[]) {
  char* filename = NULL;
  struct svm_args parameters;
  struct denseData ds;
  struct Fullproblem fp;
  struct Projected sp;

  // Input processed:
  parse_arguments(argc, argv, &filename, &parameters);
  read_file(filename, &ds);
  preprocess(&ds);

  // Projected problem size chosen temporarily
  int p = 4;


  if (parameters.test) {
    alloc_prob(&fp, ds.nInstances, p);
    init_prob(&fp, ds.nInstances, p, &ds);
    setH(&fp, &ds);
    for (int i = 0; i < fp.q; i++) {
      for (int j = 0; j < fp.p; j++) {
        std::cout << fp.partialH[i][j] << '\t';
      }
      std::cout << '\n';
    }
    for (int i = 0; i < fp.p; i++) {
      std::cout << fp.active[i] << '\n';
    }
    alloc_subprob(&sp, p, &fp, &ds);
    init_subprob(&sp, &fp, &ds);
    for (int i = 0; i < sp.p; i++) {
      std::cout << "inactive " << fp.active[i] << '\n';
      std::cout << "yaha   " << sp.yHat[i] << '\n';
      std::cout << "alpha h  "<< sp.alphaHat[i] << '\n';
      std::cout << "rHAT  "<< sp.rHat[i] << '\n';
    }
    for (int i = 0; i < sp.p; i++) {
      for (int j = i; j < sp.p; j++) {
        std::cout << sp.H[i][j] << '\t';
      }
      std::cout << '\n';
    }

    sp.alphaHat[1] = 0.5;
    sp.alphaHat[3] = 0.5;
    sp.alphaHat[2] = 4.0;
    sp.alphaHat[0] = 0.3;
    calcYTR(&sp);
    updateAlphaR(&fp, &sp);
    calculateBeta(&fp, &sp, &ds);

    int k = singleswap(&ds, &fp,&sp,2);

    for (int i = 0; i < fp.p; i++) {
      std::cout << " active of " << i << " is "<< fp.active[i] << '\n';
      std::cout << " alpha of " << i << " is "<< fp.alpha[fp.active[i]] << '\n';
    }
    return 0;


  }



  // Full problem allocated and filled in:
  alloc_prob(&fp, ds.nInstances, p);
  init_prob(&fp, ds.nInstances, p, &ds);

  // Subproblem allocated:
  alloc_subprob(&sp, p, &fp, &ds);

  // We loop until no negative entries in beta:
  int k = 1;
  int max_iters = 500;
  int itt = 0;
  int n = 0;

  while(k){

    // H matrix columns fixed and subproblem changed
    setH(&fp, &ds);
    init_subprob(&sp, &fp, &ds);

    // Conjugate gradient algorithm with interrupt if
    n = cg(&sp);

    // Computations for choosing next elements to analyze
    calcYTR(&sp);
    updateAlphaR(&fp, &sp);
    calculateBeta(&fp, &sp, &ds);
std::cout << "n is " << n << '\n';
    if (n) {
      // Case where interrupt needed:
      k = singleswap(&ds, &fp,&sp,n);
      for (int i = 0; i < sp.p; i++) {
        std::cout << "activia " << fp.active[i] << '\n';
      }
    }
    else{
      //  Without interrupt we swap in full
      k = swapMostNegative(&fp);
      for (int i = 0; i < sp.p; i++) {
        std::cout << "muller " << fp.active[i] << '\n';
      }
    }

    itt++;
    if(itt == max_iters){
      break;
    }
  }
  if (k==0) {
    std::cout << "k goes to 0." << '\n';
  }


  for (int k = 0; k < fp.n; k++) {
    std::cout << " alpha [ " << k << " ] =  " << fp.alpha[k] << '\n';
    if (fp.alpha[k]<0.0) {
      return -3;
    }
  }

  //Memory freed
  freeDenseData(&ds);
  freeFullproblem(&fp);
  freeSubProblem(&sp);
  std::cout << itt << '\n';
  return 0;
}

void   freeSubProblem( struct Projected* sp)
/* Function to free dynamically allocated memory in subproblem stuct. */
{
  free(sp->alphaHat);
  free(sp->yHat);
  free(sp->rHat);
  free(sp->H);

  free(sp->gamma);
  free(sp->rho);
  free(sp->Hrho);

  free(sp->h);
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

void  freeFullproblem(struct Fullproblem *fp)
/* Function to free dynamically allocated memory in Fullproblem struct */
{
  free(fp->alpha);
  free(fp->beta);
  free(fp->gradF);

  free(fp->active);
  free(fp->inactive);

  free(fp->partialH);
  free(fp->h);
}
