// Spatial poisson GLMM on a grid, with exponentially decaying correlation function
#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(Q);
  //DATA_MATRIX(L);
  DATA_VECTOR(wt);
  PARAMETER(ridge);
  
  PARAMETER_VECTOR(b);
  PARAMETER(ltau);
  PARAMETER_VECTOR(u);
  //PARAMETER_VECTOR(lagrange);

  using namespace density;
  int i, j, n, p;
  double r;
  
  r = asDouble(ridge);
  
  Type res=0;
  p = X.cols();
  n = X.rows();
  
  vector<Type> eta(n); 
  //eta = X*b + L*u;
  eta = X*b + u;
  
  Type logdetQ;
  
  //Eigen::Matrix<Type, Dynamic, Dynamic> diag(n-1, n-1);
  Eigen::SparseMatrix<Type> diag(n, n);
  diag.setIdentity();
  diag = ridge*diag;
  Q = Q + diag;
  
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<Type> > ldl(Q);
  //vector<Type> D=Q.ldlt().vectorD();
  vector<Type> D=ldl.vectorD();
  //for (i=0; i<n; i++)
  //  if (D[i]>0) logdetQ += log(D[i]);
  //logdetQ = log(D.head(n-1)).sum();
logdetQ = log(D).sum();
//Eigen::

  //u = u - u.mean();
  res += -Type(.5)*logdetQ + Type(0.5)*exp(ltau)*(u*(Q*u.matrix()).array()).sum() + n/2.0 * Type(log(2.0*M_PI) - ltau);
  //res += pen * abs((X.transpose()*u.matrix()).sum());
  //res += pen * abs(u.sum());
  //res += 0.1 * (u*(XXt*u.matrix()).array()).sum(); 
  //res += GMRF(Q)(u);
  //res += (u*u).sum();
  
  // logdpois = N log lam - lam
  for(i=0;i<n;i++) res -= wt[i] * (y[i]*eta[i] - exp(eta[i]));
  
  //Lagrange multipliers
  //for(j=0; j<p; j++) res -= lagrange(j) * (X.col(j).array() * u).sum() / double(n);

  return res;
}
