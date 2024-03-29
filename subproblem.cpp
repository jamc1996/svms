#include "subproblem.h"


void alloc_subprob(struct Projected *sp, int p, struct Fullproblem *fp, struct denseData *ds)
{
  sp->p = p;
  sp->C = 100.0;

  sp->alphaHat = (double*)malloc(sizeof(double)*p);
  sp->yHat = (double*)malloc(sizeof(double)*p);
  sp->rHat = (double*)malloc(sizeof(double)*p);

  sp->gamma = (double*)malloc(sizeof(double)*p);
  sp->rho = (double*)malloc(sizeof(double)*p);
  sp->Hrho = (double*)malloc(sizeof(double)*p);

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
        sp->H[i][j] += ds->data[fp->active[i]][k]*ds->data[fp->active[j]][k];
      }
      sp->H[i][j]*=ds->y[fp->active[i]]*ds->y[fp->active[j]];
    }
  }

  std::cout << sp->H[2][3] << '\n';

}

void init_symmetric(struct Projected *sp, int p)
{


  std::cout << sp->H[2][3] << '\n';

}

int cg(struct Projected *sp)
{
  double lambda, mu;
  init_error(sp);
  double rSq = inner_prod(sp->gamma,sp->gamma,sp->p);

  double newRSQ;
  int problem = 0;
  int i=0;
  int flag = 1;
  while (rSq > 0.001) {
    calc_Hrho(sp);
    for (int i = 0; i < sp->p; i++) {
  //    std::cout << "alpha is " << sp->alphaHat[i] << '\n';
    }

    for (int i = 0; i < sp->p; i++) {
  //    std::cout << "rho is " << sp->rho[i] << '\n';
    }
    std::cout << '\n';
    for (int i = 0; i < sp->p; i++) {
  //    std::cout << "Hrho is " << sp->Hrho[i] << '\n';
    }
    std::cout << '\n';
    for (int i = 0; i < sp->p; i++) {
  //    std::cout << "Gamma is " << sp->gamma[i] << '\n';
    }
std::cout << i << '\n';

    lambda = rSq/inner_prod(sp->Hrho, sp->rho, sp->p);
    std::cout << "alpah 3 is " << sp->alphaHat[3] << '\n';
    linearOp(sp->alphaHat, sp->rho, lambda, sp->p);
    std::cout << "alpah 3 is " << sp->alphaHat[3] << '\n';

    problem = checkConstraints(sp);
    if(problem){
      std::cout << "uh oh " << problem << '\n';
      linearOp(sp->alphaHat, sp->rho, -lambda, sp->p);
      std::cout << "alpah 3 is " << sp->alphaHat[3] << '\n';

      return problem;
    }

    updateGamma(sp, lambda);

    newRSQ = inner_prod(sp->gamma, sp->gamma, sp->p);

    mu = newRSQ/rSq;
    linearOp2(sp->rho, sp->gamma, mu, sp->p);

    rSq = newRSQ;
    std::cout << "Here it is " << rSq << '\n';
    i++;

  }

  for (int i = 0; i < sp->p; i++) {
    std::cout << "gamma is " << sp->gamma[i] << '\n';
  }
  return 0;
}

int checkConstraints(struct Projected* sp)
{
  double* temp = (double*)malloc(sizeof(double)*sp->p);
  constraint_projection(temp, sp->alphaHat, sp->yHat, sp->p);
  for (int i = 0; i < sp->p; i++) {
    if(temp[i]>sp->C){
      free(temp);
      return i+sp->p;
    }
    else if(temp[i]<0.0){
      free(temp);
      return i;
    }
  }
  free(temp);
  return 0;
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
void updateGamma(struct Projected *sp, double lambda)
{
  for (int i = 0; i < sp->p; i++) {
    sp->Hrho[i]*=lambda;
  }
  for (int i = 0; i < sp->p; i++) {
    sp->gamma[i] -= sp->Hrho[i];
    for (int j = 0; j < sp->p; j++) {
      sp->gamma[i] += sp->yHat[i]*sp->yHat[j]*sp->Hrho[j]/((double)sp->p);
    }
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
