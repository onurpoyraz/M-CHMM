#define ARMA_NO_DEBUG
#include <RcppArmadillo.h> 
using namespace Rcpp;
using namespace arma;
using namespace R;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List sample_latent_variable_iffbs(const field<field<Col<double>>>& pi0, 
                                  const field<field<Mat<double>>>& emission, 
                                  const field<field<field<Cube<double>>>>& beta,
                                  const Col<double>& mixing_coeficient,
                                  const Cube<int>& latent_sequences, 
                                  const Cube<int>& observation, 
                                  const int& n_group, 
                                  const int& n_states) {
    int N = observation.n_rows;
    int L = observation.n_cols;
    int C = observation.n_slices;
    int group;
    double prob;

    Mat<double> transition(n_states, n_states);
    Row<double> row_transition(n_states);
    Col<double> modifying_mass(n_states);
    Cube<int> latent(N, L, C);
    Col<double> log_likelihood(N);
    Col<double> log_lik(n_group);
    Col<int> z(N);
    Mat<double> alpha(n_states, L);
    Col<double> sampling_dist(n_states);
    Mat<int> individual(L, C);
    Row<int> covariates(C);

    auto row_softmax = [] (Row<double>& y) {
        double logsumexp = log(sum(exp(y)));
        y = exp(y - logsumexp);
    };

    auto row_logit = [] (Row<double>& y) {
        double logsumexp = log(sum(exp(y)));
        y -= logsumexp;
    };

    auto col_softmax = [] (Col<double>& y) {
        double logsumexp = log(sum(exp(y)));
        Col<double> x = exp(y - logsumexp);
        return x;
    };

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

    for(int i = 0; i < N; i++){
        individual = latent_sequences.row(i);
        // LOG LIKELIHOOD
        log_lik.zeros();
        for(int g = 0; g < n_group; g++){
            for(int chain = 0; chain < C; chain++) {
                for(int t = 0; t < L; t++){
                    if (t == 0) {
                        alpha.col(t) = pi0(g)(chain);
                    } else {
                        transition.zeros();
                        for (int c = 0; c < C; c++) {
                            if (c == chain) {
                                transition += beta(g)(chain)(chain).slice(0);
                            } else {
                                transition += beta(g)(chain)(c).slice(individual(t - 1, c));
                            }
                        }
                        transition.each_row(row_softmax);
                        alpha.col(t) = transition.t() * alpha.col(t - 1);
                    }
                    if (!IntegerVector::is_na(observation(i, t, chain))) {
                        alpha.col(t) %= emission(g)(chain).col(observation(i, t, chain));
                        prob = sum(alpha.col(t));
                        alpha.col(t) /= prob;
                        log_lik(g) += log(prob);
                    } 
                }
            }
        }
        // SAMPLE Z
        sampling_dist = softmax(log_lik);
        group = sample_latent(sampling_dist, n_group);
        z(i) = group;
        prob = dot(exp(log_lik), sampling_dist);
        log_likelihood(i) = log(prob);

        for(int chain = 0; chain < C; chain++) {
            // FORWARD_FILTERING
            for(int t = 0; t < L; t++){
                if (t == 0) {
                    alpha.col(t) = pi0(group)(chain);
                } else {
                    transition.zeros();
                    for (int c = 0; c < C; c++) {
                        if (c == chain) {
                            transition += beta(group)(chain)(chain).slice(0);
                        } else {
                            transition += beta(group)(chain)(c).slice(individual(t - 1, c));
                        }
                    }
                    transition.each_row(row_softmax);
                    alpha.col(t) = transition.t() * alpha.col(t - 1);
                }
                if (t < L - 1) {
                    covariates = individual.row(t);
                    modifying_mass = log(alpha.col(t));
                    for(int c = 0; c < C; c++) {
                        if (c != chain){
                            for(int s = 0; s < n_states; s++){
                                covariates(chain) = s;
                                transition.zeros();
                                for(int c_cov = 0; c_cov < C; c_cov++) {
                                    if (c_cov == c) {
                                        transition += beta(group)(c)(c).slice(0);
                                    } else {
                                        transition += beta(group)(c)(c_cov).slice(covariates(c_cov));
                                    }
                                }
                                row_transition = transition.row(individual(t, c));
                                row_logit(row_transition);
                                modifying_mass(s) += row_transition(individual(t+1, c));
                            }
                        }
                    }
                    alpha.col(t) = col_softmax(modifying_mass);
                }
                if (!IntegerVector::is_na(observation(i, t, chain))) {
                    alpha.col(t) %= emission(group)(chain).col(observation(i, t, chain));
                    prob = sum(alpha.col(t));
                    alpha.col(t) /= prob;
                } 
            }
            // BACKWARD SAMPLING
            for(int t = L-1; t >= 0; t--){
                if (t == L-1) {
                    sampling_dist = alpha.col(t);
                } else {
                    transition.zeros();
                    for (int c = 0; c < C; c++) {
                        if (c == chain) {
                            transition += beta(group)(chain)(chain).slice(0);
                        } else {
                            transition += beta(group)(chain)(c).slice(individual(t, c));
                        }
                    }
                    transition.each_row(row_softmax);
                    sampling_dist = alpha.col(t) % transition.col(individual(t+1, chain));
                    sampling_dist /= sum(sampling_dist);
                }
                individual(t, chain) = sample_latent(sampling_dist, n_states);
            }
        }
        latent.row(i) = individual;
    }
    return List::create(log_likelihood, latent, z);
}