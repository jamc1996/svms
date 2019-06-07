#include "subproblem.h"


void alloc_subprob(struct Projected *sp, int p, struct Fullproblem *fp, struct denseData *ds)
{
  sp->p = p;

  sp->alphaHat = (double*)malloc(sizeof(double)*p);
  sp->yHat = (double*)malloc(sizeof(double)*p);
  sp->rHat = (double*)malloc(sizeof(double)*p);

  sp->gamma = (double*)malloc(sizeof(double)*p);
  sp->rho = (double*)malloc(sizeof(double)*p);
  sp->Hrho = (double*)malloc(sizeof(double)*p);

  std::cout << "Inint" << '\n';
  //init_symmetric(sp,p);

  sp->H = (double**)malloc(sizeof(double*)*p);
  sp->h = (double*)malloc(sizeof(double)*((p*(p+1))/2));

  int j = 0;
  for (int i = 0; i < p; i++) {
    sp->H[i] = &(sp->h[j]);
    j+=(p-i-1);
  }

}

void init_subprob(struct Projected *sp, struct Fullproblem *fp, struct denseData *ds)
{
  for (int i = 0; i < sp->p; i++) {
    sp->yHat[i] = ds->y[fp->active[i]];
    sp->alphaHat[i] = fp->alpha[fp->active[i]];
    sp->rHat[i] = fp->gradF[fp->active[i]];
  }

  for (int i = 0; i < sp->p; i++) {
    for (int j = i; j < sp->p; j++) {
      sp->H[i][j] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        sp->H[i][j] += ds->data[i][k]*ds->data[j][k];
      }
    }
  }

  std::cout << sp->H[2][3] << '\n';

}

void init_symmetric(struct Projected *sp, int p)
{


  std::cout << sp->H[2][3] << '\n';

}

void cg(struct Projected *sp)
{
  double lambda, mu;
  init_error(sp);
  double rSq = inner_prod(sp->gamma,sp->gamma,sp->p);
  double newRSQ;
  while (rSq > 0.00001) {
    std::cout << "ok1" << '\n';
    calc_Hrho(sp);
    std::cout << "ok2" << '\n';

    lambda = rSq/inner_prod(sp->Hrho, sp->rho, sp->p);
    std::cout << "o32" << '\n';

    linearOp(sp->alphaHat, sp->rho, lambda, sp->p);
    std::cout << "ok4" << '\n';
    updateGamma(sp, lambda);
    std::cout << "ok5" << '\n';

    newRSQ = inner_prod(sp->gamma, sp->gamma, sp->p);
    std::cout << "ok6" << '\n';

    mu = rSq/newRSQ;
    linearOp2(sp->rho, sp->gamma, mu, sp->p);
    std::cout << "ok7" << '\n';

    break;
  }
}

void linearOp2(double* rho, double* gamma, double mu, int p)
{
  for (int i = 0; i < p; i++) {
    rho[i] *= mu;
    rho[i] += gamma[i];
  }
}

void calcYTR(struct Projected *sp)
{
  sp->ytr = 0.0;
  for (int i = 0; i < sp->p; i++) {
    sp->ytr += sp->yHat[i]*sp->rHat[i];
  }
}

void updateGamma(struct Projected *sp, double lambda)
{
  for (int i = 0; i < sp->p; i++) {
    sp->gamma[i] -= lambda*sp->Hrho[i];
    for (int j = 0; j < sp->p; j++) {
      sp->gamma[i] -= sp->yHat[i]*sp->yHat[j]*sp->Hrho[j];
    }
  }
}

void linearOp(double* alpha, double* rho, double lambda, int p)
{
  for (int i = 0; i < p; i++) {
    alpha[i] += lambda*rho[i];
  }
}

void calc_Hrho(struct Projected *sp)
{
  for (int i = 0; i < sp->p; i++) {
    sp->Hrho[i] = 0.0;
    for (int j = 0; j < i; j++) {
      sp->Hrho[i]+= sp->H[j][i]*sp->rho[j];
    }
    for (int j = i; j < sp->p; j++) {
      sp->Hrho[i]+= sp->H[i][j]*sp->rho[j];
    }
  }
}

void init_error(struct Projected* sp)
{
  constraint_projection(sp->gamma, sp->rHat, sp->yHat, sp->p);
  copy_vector(sp->rho, sp->gamma, sp->p);
}

void copy_vector(double* a, double* b, int p)
{
  for (int i = 0; i < p; i++) {
    a[i] = b[i];
  }
}

void constraint_projection(double* vecOut, double* vecIn, double* y, int p)
{
  for (int i = 0; i < p; i++) {
    vecOut[i] = vecIn[i];
    for (int j = 0; j < p; j++) {
      vecOut[i] -= vecIn[j]*(y[i]*y[j]/(double)p);
    }
  }
}

double inner_prod(double *a, double *b, int p)
{
  double val = 0.0;
  for (int i = 0; i < p; i++) {
    val+=a[i]*b[i];
  }
  return val;
}
