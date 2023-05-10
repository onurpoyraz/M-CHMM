#' A Reference Class to ChmmMcmc
#'
#' This class corresponds to CHMM training.
#' It basically manage the Chain class.
#' @field params Parameter list that defines the model parameters.
#'     MUST be supplied by user.
#'     Details of that field are in the package description.
#' @field chain_names Optional parameter to give the chains a name.
#' @field chains A list that store Chain class objects. It has the size of  number of chains.
#' @field log_likelihoods Log likelihoods of each sequence at each MCMC iteration.
#'     It is the sum of log_likelihoods of all chains.
#'     Rows corresponds to different sequences and columns corresponds to iterations
#' @field n_chains Number of chains. Class calculate it from length of observation list.
#' @field runtime Total time spend on training
#' @field diagnosis Model diagnostics and model scores.

MixtureChmm <- R6::R6Class(
  "MixtureChmm",
  portable = FALSE,
  cloneable = FALSE,
  class = FALSE,
  public = list (
    params = list(),
    diagnosis = list(),
    prior = list(),
    current = list(),
    samples = list(),
    posterior = list(),
    mcmc = list(),
    log_likelihoods = NULL,
    update_pi0 = NULL,
    update_emission = NULL,
    update_beta = NULL,
    update_latent_variable = NULL,
    update_mixing = NULL,
    test = NULL,
    initialize = function(input_observations, input_params, input_prior,
                          pi0_sampler, emission_sampler, beta_sampler, latent_variable_sampler, mixing_sampler,
                          test_input=NULL) {
      if (typeof(input_observations) == "list") {
        current$observations <<- join_matrices(input_observations)
      } else {
        current$observations <<- input_observations
      }
      params <<- input_params
      update_params()
      prior <<- input_prior
      test <<- test_input
      if (is.null(params$chain_names)) {
        params$chain_names <<- sapply(1:params$n_chains, function(x) paste("C", toString(x), sep = ""))
      }
      add_sampler("update_pi0", pi0_sampler)
      add_sampler("update_emission", emission_sampler)
      add_sampler("update_beta", beta_sampler)
      add_sampler("update_latent_variable", latent_variable_sampler)
      add_sampler("update_mixing", mixing_sampler)
    },
    add_sampler = function(name, method) {
      self[[name]] <- method
      environment(self[[name]]) <- environment(self$add_sampler)
    },
    update_params = function() {
      params$N <<- dim(current$observations)[1] ## number of sequence
      params$L <<- dim(current$observations)[2] ## length of sequence
      params$n_chains <<- dim(current$observations)[3] ## number of chains
      params$alphabet <<- as.integer(0:(params$n_states-1))
    },
    train = function() {
      ptm <- proc.time()
      for (iter in 1:params$n_iter) {
        update_mixing(iter)
        update_pi0(iter)
        update_emission(iter)
        update_beta(iter)
        update_latent_variable(iter)
      }
      finalize_training()
      diagnosis$runtime <<- proc.time() - ptm
    },
    finalize_training = function() {
      current <<- NULL
      if (params$warmup + 1 >= params$n_iter) {
        select <- 1:params$n_iter
      } else {
        select <- (params$warmup+1):params$n_iter
      }
      calculate_posterior(select)
    },
    calculate_posterior = function(select) {
      update_mixing(NULL, training = FALSE, select)
      update_pi0(NULL, training = FALSE, select)
      update_emission(NULL, training = FALSE, select)
      update_beta(NULL, training = FALSE, select)
      update_latent_variable(NULL, training = FALSE, select)
    },
    simulate = function(n_length, n_individual) {
      simulation <- simulator(posterior$pi0,
                              posterior$emission,
                              posterior$converted_beta,
                              posterior$mixing_coefficient,
                              params$n_states,
                              params$n_chains,
                              n_length,
                              n_individual)
      return (simulation)
    },
    diagnostics = function(select=(params$warmup+1):params$n_iter)
    {
      diagnosis$waic_score <<- loo::waic(t(log_likelihoods[ , select]))
      diagnosis$loo_score <<- loo::loo(t(log_likelihoods[ , select]))
    }
  )
)
