#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>   
using namespace Rcpp;
using namespace arma;
using namespace R;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Col<double> sample_mixing_coefficient(const Col<double>& prior, 
                                      const Col<int>& z) {
    Col<double> counts = prior;
    int n_individual = z.n_elem;

    auto rdirichlet = [] (Col<double>& params) {
        int n = params.n_elem;
        Col<double> sample(n);
        for(int i = 0; i < n; i++) {
            sample(i) = rgamma(params(i), 1.0);
        }
        sample /= sum(sample);
        return sample;
    };

    for(int i = 0; i < n_individual; i++) {
        counts(z(i)) += 1;
    }
    return rdirichlet(counts);
}