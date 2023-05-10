#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>   
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void cumulative_covariance(const Col<double>& data, Mat<double>& r_covariance, Col<double>& r_mean, Col<int>& r_n) {
    r_n += 1;
    Col<double> r_diff = data - r_mean;
    r_mean += r_diff / r_n(0); 
    r_covariance += ((r_n(0) - 1.) / r_n(0)) * r_diff * r_diff.t();
}