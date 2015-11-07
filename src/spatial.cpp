// Spatial poisson GLMM on a grid, with exponentially decaying correlation function
#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() ()
{
DATA_INTEGER(n);
DATA_VECTOR(y);
DATA_MATRIX(X);
Eigen::SparseMatrix<Type>(Q);
DATA_VECTOR(wt);
PARAMETER_VECTOR(b);
PARAMETER(log_sigma);
PARAMETER_VECTOR(u);
Type sigma2=exp(2.0*log_sigma);

using namespace density;
int i,j;
Type res=0;

vector<Type> eta(n); 
eta = X*b + exp(log_sigma)*u;

GMRF_t<Type> neg_log_density(Q);
res += neg_log_density(u);

// logdpois = N log lam - lam
for(i=0;i<n;i++) res -= wt[i] * (y[i]*eta[i]-exp(eta[i]));

return res;

  }
