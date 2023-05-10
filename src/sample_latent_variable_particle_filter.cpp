#define ARMA_NO_DEBUG
#include <RcppArmadillo.h> 
using namespace Rcpp;
using namespace arma;
using namespace R;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List sample_latent_variable_particle_filter(const field<field<Col<double>>>& pi0,
                                            const field<field<Mat<double>>>& emission, 
                                            const field<field<field<Cube<double>>>>& beta,
                                            const Col<double>& mixing_coeficient,  
                                            const Cube<int>& observation, 
                                            const int& n_group,
                                            const int& n_states, 
                                            const int& P) {
    int N = observation.n_rows;
    int L = observation.n_cols;
    int C = observation.n_slices;
    int sample, trajectory, previous_t, current_state, group;
    double prob, ess;
    Col<uword> resampling(P);
    
    Row<double> row_transition(n_states);
    Col<double> log_likelihood(N);
    Col<double> log_lik(n_group);
    Col<int> z(N);
    Cube<int> latent(N, L, C);
    Cube<int> individual_particles(L, C, P);
    Cube<int> group_particles(L, C, n_group);
    Mat<double> weights(P, L, fill::zeros);
    Mat<double> weights_on_time(P, C);
    Col<double> weight_sum(P);
    Col<double> weight_prob(P);
    Col<double> sampling_dist(n_states);

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

    auto softmax_z = [mixing_coeficient] (Col<double>& x) {
        Col<double> log_prob = x + mixing_coeficient;
        double logsumexp = log(sum(exp(log_prob)));
        Col<double> prob = exp(log_prob - logsumexp);
        return prob;
    };

    auto softmax = [] (Col<double>& log_prob) {
        double logsumexp = log(sum(exp(log_prob)));
        Col<double> prob = exp(log_prob - logsumexp);
        return prob;
    };

    auto row_softmax = [] (Row<double>& y) {
        double logsumexp = log(sum(exp(y)));
        y = exp(y - logsumexp);
    };

    auto systematic_resampling = [] (Col<double>& weight_prob, const int& size) {
        Col<double> cdf = cumsum(weight_prob);
        Col<double> uvector(size);
        Col<uword> out(size);
        uword intermediate = 0;
        double interval = 1.0 / size;
        for (uword i = 0; i < size; i++){
            if (i == 0) {
                uvector(i) = runif(0.0, interval);
            } else {
                uvector(i) = uvector(i-1) + interval;
            }
            for (uword s = intermediate; s < size; s++) {    
                if (uvector(i) <= cdf(s)) {
                    out(i) = s;
                    intermediate = s;
                    break;
                }
            }
        }
        return out;
    };

    for(int i = 0; i < N; i++){
        log_lik.zeros();
        // FORWARD_SAMPLING
        for(int g = 0; g < n_group; g++){
            previous_t = 0;
            for(int t = 0; t < L; t++){
                // RESAMPLING PART
                if (t > 1) {
                    weight_sum = sum(weights.cols(previous_t, t-1), 1);
                    weight_prob = softmax(weight_sum);
                    ess = 1.0 / (P * accu(pow(weight_prob, 2)));
                    if (ess < 0.5) {
                        previous_t = t;
                        resampling = systematic_resampling(weight_prob, P);
                        individual_particles = individual_particles.slices(resampling);
                    }
                }
                for (int particle = 0; particle < P; particle++) {
                    for(int chain = 0; chain < C; chain++) {
                        if (t == 0) {
                            sampling_dist = pi0(g)(chain);
                        } else {
                            row_transition.zeros();
                            current_state = individual_particles(t - 1, chain, particle);
                            for (int c = 0; c < C; c++) {
                                if (c == chain) {
                                    row_transition += beta(g)(chain)(chain).slice(0).row(current_state);
                                } else {
                                    row_transition += beta(g)(chain)(c).slice(individual_particles(t-1, c, particle)).row(current_state);
                                }
                            }
                            row_softmax(row_transition);
                            sampling_dist = row_transition.t();
                        }
                        if (!IntegerVector::is_na(observation(i, t, chain))) {
                            sampling_dist %= emission(g)(chain).col(observation(i, t, chain));
                            prob = sum(sampling_dist);
                            sampling_dist /= prob;
                        } else {
                            prob = 1;
                        }
                        weights_on_time(particle, chain) = log(prob);
                        sample = sample_latent(sampling_dist, n_states);
                        individual_particles(t, chain, particle) = sample;
                    }
                }
                weights.col(t) = sum(weights_on_time, 1);
            }
            log_lik(g) = sum(log(mean(exp(weights), 0)));
            weight_sum = sum(weights.cols(previous_t, L-1), 1);
            weight_prob = softmax(weight_sum);
            trajectory = sample_latent(weight_prob, P);
            group_particles.slice(g) = individual_particles.slice(trajectory);
        }

        // SAMPLE Z
        sampling_dist = softmax_z(log_lik);
        group = sample_latent(sampling_dist, n_group);
        z(i) = group;
        prob = dot(exp(log_lik), sampling_dist);
        log_likelihood(i) = log(prob);

        latent.row(i) = group_particles.slice(group);
    }
    return List::create(log_likelihood, latent, z);
}

