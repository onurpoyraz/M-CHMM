#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>   
using namespace Rcpp;
using namespace arma;
using namespace R;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Mat<double> sample_emission(const Mat<double> prior, 
                            const Cube<int>& latent_sequences, 
                            const Cube<int>& observations, 
                            const Col<int>& z,
                            const int& n_states, 
                            const int& n_observation, 
                            const int& group,
                            const int& chain) {
    int x, y;
    int n_row = latent_sequences.n_rows;
    int n_col = latent_sequences.n_cols;
    Mat<double> counts(n_states, n_observation);

    auto rdirichlet = [] (Mat<double>& params) {
        int n_row = params.n_rows;
        int n_col = params.n_cols;
        Mat<double> sample(n_row, n_col);
        for(int i = 0; i < n_row; i++) {
            for(int j = 0; j < n_col; j++) {
                sample(i, j) = rgamma(params(i, j), 1.0);
            }
            sample.row(i) /= sum(sample.row(i));
        }
        return sample;
    };

    counts = prior;
    for(int i = 0; i < n_row; i++) {
        if (z(i) == group) {
            for(int j = 0; j < n_col; j++) {
                y = observations(i, j, chain);
                if(!IntegerVector::is_na(y)) {
                    x = latent_sequences(i, j, chain);
                    counts(x, y) += 1.0;
                }
            }
        }
    }
    return rdirichlet(counts);
}
