#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Col<double> multivariate_normal(const Col<double>& mean, const Mat<double>& covariance) {
    Col<double> sample = mvnrnd(mean, covariance, 1);
    return sample;
}
