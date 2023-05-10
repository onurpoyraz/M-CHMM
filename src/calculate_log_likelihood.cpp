#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>   
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double calculate_log_likelihood(const Cube<int>& latent, 
                                const field<Cube<double>>& beta, 
                                const Col<int>& z,
                                const int& n_states,
                                const int& group,
                                const int& chain) {
    int N = latent.n_rows;
    int L = latent.n_cols;
    int C = latent.n_slices;
    double log_prob = 0.0;
    Mat<int> individual (L, C);
    Mat<double> log_transition (n_states, n_states);
    Row<double> row_transition (n_states);

    auto row_logit = [] (Row<double>& y) {
        double logsumexp = log(sum(exp(y)));
        y -= logsumexp;
    };

    for(int i = 0; i < N; i++) {
        if (z(i) == group) {
            individual = latent.row(i);
            for(int t = 1; t < L; t++) {
                log_transition.zeros();
                for (int c = 0; c < C; c++) {
                    if (c == chain) {
                        log_transition += beta(chain).slice(0);
                    } else {
                        log_transition += beta(c).slice(individual(t - 1, c));
                    }
                }
                row_transition = log_transition.row(individual(t - 1, chain));
                row_logit(row_transition);
                log_prob += row_transition(individual(t, chain));
            }
        }
    }
    return log_prob;
}