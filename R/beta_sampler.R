beta_sampler <- function(iter, training = TRUE, select = NULL) {
  calculate_log_posterior <- function(converted_beta_values, beta_values, group, chain) {
    log_likelihood <- calculate_log_likelihood(current$latent_sequences, converted_beta_values, current$z, params$n_states, as.integer(group - 1), as.integer(chain - 1))
    log_prior <- calculate_log_prior(beta_values)
    log_posterior <- log_likelihood + log_prior
    return(log_posterior)
  }
  calculate_log_prior <- function(beta_values) {
    log_prior_intercept <- sum(dnorm(beta_values[1:params$dim_intercept], 0, params$sd_intercept, log = TRUE))
    log_prior_coefficient <- prior_density(beta_values[-(1:params$dim_intercept)], params$prior_distribution, params$sd_coefficient)
    log_prior <- log_prior_intercept + log_prior_coefficient
    return(log_prior)
  }
  propose <- function (iter, group, chain) {
    proposal <- list()
    if (iter <= params$warmup) {
      if (iter > as.integer(params$warmup/2)) {
        cumulative_covariance(current$beta[[group]][[chain]], mcmc$covariance[[group]][[chain]], mcmc$mean[[group]][[chain]], mcmc$n[[group]][[chain]])
      }
      covariance <- diag(params$dim_beta) * (mcmc$c[[group]][[chain]]^2 * params$step_size[[chain]])
    } else {
      cumulative_covariance(current$beta[[group]][[chain]], mcmc$covariance[[group]][[chain]], mcmc$mean[[group]][[chain]], mcmc$n[[group]][[chain]])
      covariance <- mcmc$covariance[[group]][[chain]] * (mcmc$c[[group]][[chain]]^2 / (mcmc$n[[group]][[chain]] - 1))
    }
    proposal$beta <- multivariate_normal(current$beta[[group]][[chain]], covariance)
    if (iter < params$covariate_lag) {
      proposal$beta[-(1:params$dim_intercept)] <- as.numeric(0)
    }
    proposal$converted_beta <- convert_beta(proposal$beta, params$n_states, params$n_chains, as.integer(chain - 1))
    return(proposal)
  }
  metropolis <- function (proposal, group, chain) {
    previous_likelihood <- calculate_log_posterior(current$converted_beta[[group]][[chain]],
                                                   current$beta[[group]][[chain]],
                                                   group,
                                                   chain)
    proposal_likelihood <- calculate_log_posterior(proposal$converted_beta,
                                                   proposal$beta,
                                                   group,
                                                   chain)
    ratio <- exp(proposal_likelihood - previous_likelihood)
    acceptance_rate <- runif(1)
    if (ratio > acceptance_rate) {
      return (TRUE)
    } else {
      return (FALSE)
    }
  }
  update_samples <- function(iter) {
    for (group in 1:params$n_group) {
      if (iter == 1) {
        samples$beta[[group]] <<- list()
        current$beta[[group]] <<- list()
        current$converted_beta[[group]] <<- list()
        mcmc$covariance[[group]] <<- list()
        mcmc$mean[[group]] <<- list()
        mcmc$n[[group]] <<- list()
        mcmc$beta_acceptance[[group]] <<- list()
        mcmc$beta_acceptance_warmup[[group]] <<- list()
        mcmc$c[[group]] <<- list()
        mcmc$adaptation[[group]] <<- list()
        params$dim_intercept <<- params$n_states * (params$n_states - 1)
        params$dim_coefficient <<- params$n_covariates * params$n_states * (params$n_states - 1)^2
        params$dim_beta <<- params$dim_intercept + params$dim_coefficient
        for (chain in 1:params$n_chains) {
          samples$beta[[group]][[chain]] <<- array(as.numeric(NA), c(params$dim_beta, params$n_iter))
          current$beta[[group]][[chain]] <<- array(as.numeric(0), c(params$dim_beta))
          current$converted_beta[[group]][[chain]] <<- convert_beta(current$beta[[group]][[chain]], params$n_states, params$n_chains, as.integer(chain - 1))
          samples$beta[[group]][[chain]][ , iter] <<- current$beta[[group]][[chain]]
          mcmc$covariance[[group]][[chain]] <<- array(as.numeric(0), c(params$dim_beta, params$dim_beta))
          mcmc$mean[[group]][[chain]] <<- array(as.numeric(0), c(params$dim_beta))
          mcmc$n[[group]][[chain]] <<- as.integer(0)
          mcmc$beta_acceptance[[group]][[chain]] <<- as.numeric(0)
          mcmc$beta_acceptance_warmup[[group]][[chain]] <<- as.numeric(0)
          mcmc$c[[group]][[chain]] <<- 2.38 / sqrt(params$dim_beta)
          mcmc$adaptation[[group]][[chain]] <<- mcmc$c[[group]][[chain]] / 100
        }
      } else {
        for (chain in 1:params$n_chains) {
          if (iter == (params$warmup + 1)) {
            mcmc$c[[group]][[chain]] <<- 2.38 / sqrt(params$dim_beta)
          }
          proposal <- propose(iter, group, chain)
          acceptance <- metropolis(proposal, group, chain)
          if (acceptance) {
            if (iter > params$warmup) {
              mcmc$c[[group]][[chain]] <<- mcmc$c[[group]][[chain]] + 2.3 * (mcmc$adaptation[[group]][[chain]] / sqrt(iter - params$warmup))
              mcmc$beta_acceptance[[group]][[chain]] <<- mcmc$beta_acceptance[[group]][[chain]] + 1
            } else {
              mcmc$c[[group]][[chain]] <<- mcmc$c[[group]][[chain]] + 2.3 * (mcmc$adaptation[[group]][[chain]] / sqrt(iter))
              mcmc$beta_acceptance_warmup[[group]][[chain]] <<- mcmc$beta_acceptance_warmup[[group]][[chain]] + 1
            }
            current$beta[[group]][[chain]] <<- proposal$beta
            current$converted_beta[[group]][[chain]] <<- proposal$converted_beta
          } else {
            if (iter > params$warmup) {
              mcmc$c[[group]][[chain]] <<- mcmc$c[[group]][[chain]] - (mcmc$adaptation[[group]][[chain]] / sqrt(iter - params$warmup))
            } else {
              mcmc$c[[group]][[chain]] <<- mcmc$c[[group]][[chain]] - (mcmc$adaptation[[group]][[chain]] / sqrt(iter))
            }
          }
          samples$beta[[group]][[chain]][ , iter] <<- current$beta[[group]][[chain]]
      }
    }
    }
  }
  update_posterior <- function(select) {
    for (group in 1:params$n_group) {
      posterior$beta[[group]] <<- list()
      posterior$converted_beta[[group]] <<- list()
      for (chain in 1:params$n_chains) {
        posterior$beta[[group]][[chain]] <<- rowMeans(samples$beta[[group]][[chain]][ , select])
        posterior$converted_beta[[group]][[chain]] <<- convert_beta(posterior$beta[[group]][[chain]], params$n_states, params$n_chains, as.integer(chain - 1))
        mcmc$beta_acceptance[[group]][[chain]] <<- mcmc$beta_acceptance[[group]][[chain]] / (params$n_iter - params$warmup)
        mcmc$beta_acceptance_warmup[[group]][[chain]] <<- mcmc$beta_acceptance_warmup[[group]][[chain]] / params$warmup
      }
      names(samples$beta[[group]]) <<- params$chain_names
      names(posterior$beta[[group]]) <<- params$chain_names
      names(posterior$converted_beta[[group]]) <<- params$chain_names
    }
  }
  if (training) {
    update_samples(iter)
  } else {
    update_posterior(select)
  }
}
