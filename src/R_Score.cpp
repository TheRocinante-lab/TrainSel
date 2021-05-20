#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Dense>


// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::DiagonalMatrix;


// [[Rcpp::export]]
float r_score(Eigen::MatrixXd x, Eigen::MatrixXd x0) {
  int nr = x.rows();
  int nc = x.cols();
  int n0r = x0.rows();
  Eigen::MatrixXd A  = x.transpose()*((x*x.transpose()+MatrixXd::Identity(nr, nr)*(1.0/nc)).inverse());
  Eigen::MatrixXd IJ = MatrixXd::Identity(n0r, n0r)-(MatrixXd::Identity(n0r, n0r)*(1.0/n0r));
  float q1 = n0r-1+(IJ*x0).array().square().sum();
  Eigen::MatrixXd IJX0A = IJ*x0*A;
  Eigen::MatrixXd IJX0AX = IJX0A*x;
  float q2 = IJX0A.array().square().sum() + IJX0AX.array().square().sum();
  float q12 = (x0.transpose()*IJX0AX).trace();
  return q12/(sqrt(q1*q2));
}
