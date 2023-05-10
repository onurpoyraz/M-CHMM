#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Cube<int> join_matrices(const field<Mat<int>>& x) {
    int n = x.n_elem;
    int n_row = x(0).n_rows;
    int n_col = x(0).n_cols;
    Cube<int> output (n_row, n_col, n);
    for(int i = 0; i < n; i++) {
        output.slice(i) = x(i);
    }
    return output;
}