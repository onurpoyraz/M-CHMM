#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>   
using namespace Rcpp;
using namespace arma;
using namespace R;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Col<double> sample_pi0(const Col<double>& prior, 
                       const Cube<int>& latent_sequences,
                       const Col<int>& z,
                       const int& n_states,  
                       const int& group,
                       const int& chain) {
    uvec selector = find(z == group);
    Col<int> initials = latent_sequences.slice(chain).col(0);
    Col<int> idx = initials(selector);
    Col<double> counts(n_states);
    int n_individual = idx.n_elem;

    auto rdirichlet = [] (Col<double>& params) {
        int n = params.n_elem;
        Col<double> sample(n);
        for(int i = 0; i < n; i++) {
            sample(i) = rgamma(params(i), 1.0);
        }
        sample /= sum(sample);
        return sample;
    };

    counts = prior;
    for(int i = 0; i < n_individual; i++) {
        counts(idx(i)) += 1;
    }
    return rdirichlet(counts);
}