#define ARMA_NO_DEBUG
#include <RcppArmadillo.h> 
using namespace Rcpp;
using namespace arma;
using namespace R;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List sample_latent_variable_vsffbs(const field<field<Col<double>>>& pi0, 
                                   const field<field<Mat<double>>>& emission, 
                                   const field<field<field<Cube<double>>>>& beta,
                                   const Col<double>& mixing_coeficient,  
                                   const Cube<int>& observation, 
                                   const int& n_group,
                                   const int& n_states) {
    int N = observation.n_rows;
    int L = observation.n_cols;
    int C = observation.n_slices;
    int group;
    double prob;
    
    Col<double> log_likelihood(N);
    Col<double> log_lik(n_group);
    Col<int> z(N);
    Cube<int> latent(N, L, C);
    Cube<double> alpha(n_states, L, C);
    Mat<double> transition(n_states, n_states);
    Cube<double> effect(n_states, n_states, C);
    field<Cube<double>> alpha_field(n_group);
    Col<double> back_transition(n_states);
    Col<double> sampling_dist(n_states);
    Row<double> row_transition(n_states);

    auto sample_latent = [] (Col<double>& sampling_distribution, const int& n_states) {
        Col<double> cdf = cumsum(sampling_distribution);
        double u = runif(0.0, 1.0);
        int s;
        for (s = 0; s < n_states - 1; s++) {    
            if (u <= cdf(s)) {
                return s;
            }
        }
        return s;
    };

    auto softmax = [mixing_coeficient] (Col<double>& x) {
        Col<double> log_prob = x + mixing_coeficient;
        double logsumexp = log(sum(exp(log_prob)));
        Col<double> prob = exp(log_prob - logsumexp);
        return prob;
    };

    auto row_softmax = [] (Row<double>& y) {
        double logsumexp = log(sum(exp(y)));
        y = exp(y - logsumexp);
    };
    auto col_softmax = [] (Col<double>& y) {
        double logsumexp = log(sum(exp(y)));
        y = exp(y - logsumexp);
    };
    auto row_logit = [] (Row<double>& y) {
        double logsumexp = log(sum(exp(y)));
        y -= logsumexp;
    };

    for(int i = 0; i < N; i++){
        log_lik.zeros();
        for(int g = 0; g < n_group; g++){
            for(int t = 0; t < L; t++){
                // FORWARD_FILTERING
                for(int chain = 0; chain < C; chain++) {
                    if (t == 0) {
                        alpha.slice(chain).col(t) = pi0(g)(chain);
                    } else {
                        effect.zeros();
                        for(int c = 0; c < C; c++) {
                            if (c != chain){
                                for(int s = 0; s < n_states; s++){
                                    effect.slice(c) += beta(g)(chain)(c).slice(s) * alpha(s , t-1, c);
                                }
                            }
                            else{
                                effect.slice(c) = beta(g)(chain)(chain).slice(0);
                            }
                        }
                        transition = sum(effect, 2);
                        transition.each_row(row_softmax);
                        alpha.slice(chain).col(t) = transition.t() * alpha.slice(chain).col(t - 1);
                    }
                    if (!IntegerVector::is_na(observation(i, t, chain))) {
                        alpha.slice(chain).col(t) %= emission(g)(chain).col(observation(i, t, chain));
                        prob = sum(alpha.slice(chain).col(t));
                        alpha.slice(chain).col(t) /= prob;
                        log_lik(g) += log(prob);
                    } 
                }
            }
            alpha_field(g) = alpha;
        }
        // SAMPLE Z
        sampling_dist = softmax(log_lik);
        group = sample_latent(sampling_dist, n_group);
        z(i) = group;
        prob = dot(exp(log_lik), sampling_dist);
        log_likelihood(i) = log(prob);
        alpha = alpha_field(group);
        for(int t = L-1; t >= 0; t--){
            // BACKWARD SAMPLING
            for(int chain = 0; chain < C; chain++) {
                if (t == L-1) {
                    sampling_dist = alpha.slice(chain).col(t);
                } else {
                    back_transition.zeros();
                    for(int c = 0; c < C; c++) {
                        effect.zeros();
                        for(int c_cov = 0; c_cov < C; c_cov++) {
                            if (c_cov != c){
                                for(int s = 0; s < n_states; s++){
                                    effect.slice(c_cov) += beta(group)(c)(c_cov).slice(s) * alpha(s , t, c_cov);
                                }
                            }
                            else {
                                effect.slice(c_cov) = beta(group)(c)(c).slice(0);
                            }
                        }
                        if (c == chain) {
                            transition = sum(effect, 2);
                            transition.each_row(row_logit);
                            back_transition += transition.col(latent(i, t+1, chain));
                        }
                        else {
                            for(int s = 0; s < n_states; s++){
                                effect.slice(chain) = beta(group)(c)(chain).slice(s);
                                transition = sum(effect, 2);
                                transition.each_row(row_logit);
                                back_transition(s) += log(sum(exp(transition.col(latent(i, t+1, c))) % alpha.slice(c).col(t)));
                            }
                        }
                    }
                    col_softmax(back_transition);
                    sampling_dist = alpha.slice(chain).col(t) % back_transition;
                    sampling_dist /= sum(sampling_dist);
                }
                latent(i, t, chain) = sample_latent(sampling_dist, n_states);
            }
        }
    }
    return List::create(log_likelihood, latent, z);
}

