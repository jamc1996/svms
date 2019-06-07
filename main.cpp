#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"

#include <iostream>

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
  std::cout << "ok" << '\n';


  alloc_subprob(&subP, p, &fullP, &ds);
  init_subprob(&subP, &fullP, &ds);

  cg(&subP);
  for (int i = 0; i < p; i++) {
    std::cout << subP.alphaHat[i] << '\n';
  }


  updateAlphaR(&fullP, &subP);


  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
  //    std::cout << subP.H[i][j] << '\t';
    }
    std::cout << '\n';
  }



  return 0;
}
