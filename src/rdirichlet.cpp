#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace R;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Mat<double> mat_rdirichlet(Mat<double>& x) {
    int n_row = x.n_rows;
    int n_col = x.n_cols;
    Mat<double> sample(n_row, n_col);
    for(int i = 0; i < n_row; i++) {
        for(int j = 0; j < n_col; j++) {
            sample(i, j) = rgamma(x(i, j), 1.0);
        }
        sample.row(i) /= sum(sample.row(i));
    }
    return sample;
}


// [[Rcpp::export]]
Col<double> vec_rdirichlet(Col<double>& x) {
    int n = x.n_elem;
    Col<double> sample(n);
    for(int i = 0; i < n; i++) {
        sample(i) = rgamma(x(i), 1.0);
    }
    sample /= sum(sample);
    return sample;
}