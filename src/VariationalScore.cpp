#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
NumericVector VariationalVar(const Eigen::MatrixXd cholV, Eigen::MatrixXd S)
{
  int n = S.rows();
  NumericVector v(n);
  
  for (int i=0; i<n; i++)
    v(i) = (S.row(i) * cholV).array().pow(2).sum();
  
  return v;
}


// [[Rcpp::export]]
NumericVector VariationalVarIndep(const Eigen::VectorXd diagV, Eigen::MatrixXd S)
{
  int n = S.rows();
  NumericVector v(n);
  
  for (int i=0; i<n; i++)
    v(i) = (S.row(i).array().pow(2) * diagV.array()).sum();
  
  return v;
}




// [[Rcpp::export]]
Eigen::MatrixXd VariationalScore(Eigen::VectorXd mu, Eigen::VectorXd wt, double tau, Eigen::VectorXd v, const Eigen::MatrixXd V, const Eigen::MatrixXd S)
{
  int p = S.cols();
  
  // Gradient of the lower likelihood bound:
  Eigen::VectorXd diag = mu.array() * wt.array() * v.array();
  Eigen::MatrixXd grad = -S.transpose() * diag.asDiagonal() * S;
  
  grad += V.inverse();
  for (int i=0; i<p; i++)
    grad.diagonal()[i] -= tau;
  grad *= 0.5;
  
  return grad;
}
