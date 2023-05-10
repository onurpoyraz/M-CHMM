#define ARMA_NO_DEBUG
#include <RcppArmadillo.h> 
using namespace Rcpp;
using namespace arma;
using namespace R;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
List simulator(const field<field<Col<double>>>& pi0, 
               const field<field<Mat<double>>>& emission, 
               const field<field<field<Cube<double>>>>& beta, 
               const Col<double>& mixing,
               const int& n_states, 
               const int& n_chains,
               const int& n_length, 
               const int& n_individual) {
    int sample, current_state, group;
    Row<double> row_transition(n_states);
    Cube<int> latent(n_individual, n_length, n_chains);
    Cube<int> observation(n_individual, n_length, n_chains);
    Col<int> z(n_individual);
    Mat<int> individual_latent(n_length, n_chains);
    Mat<int> individual_observation(n_length, n_chains);
    Mat<int> joint_latent(n_individual, n_length, fill::ones); // to be deleted
    Mat<int> joint_observation(n_individual, n_length, fill::ones); // to be deleted

    auto draw_sample = [] (const Col<double>& sampling_distribution) {
        Col<double> cdf = cumsum(sampling_distribution);
        int n_states = cdf.n_elem;
        double u = runif(0.0, 1.0);
        int s;
        for (s = 0; s < n_states - 1; s++) {    
            if (u <= cdf(s)) {
                return s;
            }
        }
        return s;
    };

    auto row_softmax = [] (Row<double>& y) {
        double logsumexp = log(sum(exp(y)));
        y = exp(y - logsumexp);
    };

    for(int i = 0; i < n_individual; i++){
        group = draw_sample(mixing);
        z(i) = group;
        for(int t = 0; t < n_length; t++){
            for(int chain = 0; chain < n_chains; chain++) {
                if (t == 0) {
                    sample = draw_sample(pi0(group)(chain));
                    individual_latent(t, chain) = sample;
                }
                else {
                    row_transition.zeros();
                    current_state = individual_latent(t - 1, chain);
                    for (int c = 0; c < n_chains; c++) {
                        if (c == chain) {
                            row_transition += beta(group)(chain)(chain).slice(0).row(current_state);
                        } else {
                            row_transition += beta(group)(chain)(c).slice(individual_latent(t - 1, c)).row(current_state);
                        }
                    }
                    row_softmax(row_transition);
                    sample = draw_sample(row_transition.t());
                    individual_latent(t, chain) = sample;
                }
                individual_observation(t, chain) = draw_sample(emission(group)(chain).row(sample).t());
            }
            if (all(individual_latent.row(t) == 0)){ // to be deleted
                joint_latent(i, t) = 0; // to be deleted
            } // to be deleted
            if (all(individual_observation.row(t) == 0)){ // to be deleted
                joint_observation(i, t) = 0; // to be deleted
            } // to be deleted
        }
        latent.row(i) = individual_latent;
        observation.row(i) = individual_observation;
    }
    List output;
    output["z"] = z;
    output["latent_sequences"] = latent;
    output["observations"] = observation;
    output["joint_latent_sequences"] = joint_latent; // to be deleted
    output["joint_observations"] = joint_observation; // to be deleted
    return output;
}