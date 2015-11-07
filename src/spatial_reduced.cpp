// Spatial poisson GLMM on a grid, with exponentially decaying correlation function
#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(Q);
  DATA_SCALAR(logdetQ);
  DATA_MATRIX(SRE);
  DATA_VECTOR(wt);
  DATA_SCALAR(ltau);

  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(u);

  using namespace density;
  using namespace Eigen;
  
  int i, j, n, p, r;
  
  Type res=0;
  p = X.cols();
  n = X.rows();
  r = SRE.ncol();
  
  vector<Type> eta(n); 
  eta = X*b + SRE*u;
  
  res += -Type(.5)*logdetQ + Type(0.5)*exp(ltau)*(u.array()*(Q*u.matrix()).array()).sum() + r/2.0 * Type(log(2.0*M_PI) - ltau);

  for(i=0;i<n;i++) res -= wt[i] * (y[i]*eta[i] - exp(eta[i]));
  
  return res;
}
