#define ARMA_NO_DEBUG
#include <RcppArmadillo.h> 
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List convert_beta(const Col<double>& beta_values, 
                  const int& n_states, 
                  const int& n_chain,
                  const int& chain) {
    List beta(n_chain);
    Cube<double> beta_cube(n_states, n_states, n_states);
    Mat<double> beta_intercept(n_states, n_states);
    int iterator = 0;

    beta_intercept.zeros();
    for(int j = 0; j < n_states - 1; j++){
        for(int i = 0; i < n_states; i++){
            beta_intercept(i, j) = beta_values(iterator);
            iterator++;
        }
    }
    beta_intercept.col(n_states - 1) = -sum(beta_intercept, 1);

    beta_cube.zeros();
    for(int k = 0; k < n_states; k++) {
        beta_cube.slice(k) = beta_intercept;
    }
    beta[chain] = beta_cube;
    
    for(int c = 0; c < n_chain; c++) {
        if (c != chain){
            beta_cube.zeros();
            for(int k = 1; k < n_states; k++) {
                for(int j = 0; j < n_states - 1; j++){
                    for(int i = 0; i < n_states; i++){
                        beta_cube(i, j, k) = beta_values(iterator);
                        iterator++;
                    }
                }
                beta_cube.slice(k).col(n_states - 1) = -sum(beta_cube.slice(k), 1);
            }
            beta[c] = beta_cube;
        }
    }
    return beta;
}