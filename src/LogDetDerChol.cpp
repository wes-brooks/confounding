#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
NumericVector LogDetDerChol(const Eigen::MatrixXd L, const Eigen::MatrixXd S, const Eigen::MatrixXd X, const Eigen::VectorXd mu, const Eigen::VectorXd wt) 
{
  int r = S.cols();
  int p = X.cols();
  
  // dM/d(beta): the change in the Hessian matrix with change in beta.
  Eigen::MatrixXd *gr = new Eigen::MatrixXd[p];
  for (int i=0; i<p; i++)
    gr[i] = Eigen::MatrixXd::Zero(r, r);
  
  Eigen::MatrixXd StS = Eigen::MatrixXd::Zero(r, r);
  for (int j=0; j<p; j++)
  {
    VectorXd Xcol = X.col(j);
    Eigen::VectorXd diag = wt.array() * mu.array() * Xcol.array();
    gr[j] = S.transpose() * diag.asDiagonal() * S;
  }
  
  NumericVector dL(p+1);
  Eigen::MatrixXd F;
  for (int l=0; l<r; l++)
  {
    //if (verbose) cat('.') 
    
    //-----------------------------------------
    // Calculate the derivative of the log-determinant of the Hessian's Cholesky factorization, following reverse-differentiation steps on page 138 of
    // S.P. Smith paper (Journal of Computational and Graphical Statistics, 1995, 4(2)):
    F = Eigen::MatrixXd::Zero(r, r);
    F(l,l) = 1;
    int lim = min(r-2, l);
    
    for(int k=lim; k>=0; k--) {
      for(int j=k+1; j<=l+1 && j<r; j++) {
        for(int i=j; i<=l && i<r; i++) {
          F(i,k) -= F(i,j) * L(j,k);
          F(j,k) -= F(i,j) * L(i,k);
        }
        F(j,k) = F(j,k) / L(k, k);
        F(k,k) -= L(j,k) * F(j,k);
      }
      F(k,k) = F(k,k) / L(k,k) / 2;
    }
    
    dL(p) += F.diagonal().sum() / L(l,l);
    for (int j=0; j<p; j++)
      dL(j) += (F.array()*gr[j].array()).sum();
  }

  
  return(dL);
}


// [[Rcpp::export]]
NumericVector LogDetDerChol2(const Eigen::MatrixXd L, const Eigen::MatrixXd S, const Eigen::MatrixXd X, const Eigen::VectorXd mu, const Eigen::VectorXd wt, double tau) 
{
  int r = S.cols();
  int p = X.cols();
  
  // dM/d(beta): the change in the Hessian matrix with change in beta.
  Eigen::MatrixXd *gr = new Eigen::MatrixXd[p];
  for (int i=0; i<p; i++)
    gr[i] = Eigen::MatrixXd::Zero(r, r);
  
  Eigen::MatrixXd StS = Eigen::MatrixXd::Zero(r, r);
  for (int j=0; j<p; j++)
  {
    VectorXd Xcol = X.col(j);
    Eigen::VectorXd diag = wt.array() * mu.array() * Xcol.array();
    gr[j] = S.transpose() * diag.asDiagonal() * S;
  }
  
  NumericVector dL(p+1);
  Eigen::MatrixXd F;
  F = Eigen::MatrixXd::Zero(r, r);
  for (int l=0; l<r; l++)
    F(l,l) = 1/L(l,l);
  

  //if (verbose) cat('.') 
  
  //-----------------------------------------
  // Calculate the derivative of the log-determinant of the Hessian's Cholesky factorization, following reverse-differentiation steps on page 138 of
  // S.P. Smith paper (Journal of Computational and Graphical Statistics, 1995, 4(2)):
  
  int lim = r-2;
  
  for(int k=r-2; k>=0; k--) {
    for(int j=k+1; j<r; j++) {
      for(int i=j; i<r; i++) {
        F(i,k) -= F(i,j) * L(j,k);
        F(j,k) -= F(i,j) * L(i,k);
      }
      F(j,k) = F(j,k) / L(k, k);
      F(k,k) -= L(j,k) * F(j,k);
    }
    F(k,k) = F(k,k) / L(k,k) / 2;
  }
  
  dL(p) = tau * F.diagonal().sum();
  for (int j=0; j<p; j++)
    dL(j) = (F.array()*gr[j].array()).sum();
  
  return(dL);
}