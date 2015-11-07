#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
VectorXd ldl(const Eigen::SparseMatrix<double> Q) 
{
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > ldl(Q);
  return(ldl.vectorD());
}