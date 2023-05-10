latent_variable_vsffbs_sampler <- function (iter, training = TRUE, select = NULL) {
  update_samples <- function(iter) {
    if (iter == 1) {
      samples$z <<- array(as.integer(NA), c(params$N, params$n_iter))
      log_likelihoods <<- array(as.numeric(NA), c(params$N, params$n_iter))
      posterior$latent_sequences <<- array(as.integer(0), c(params$N, params$L, params$n_chains))
    }
    output <- sample_latent_variable_vsffbs(current$pi0,
                                            current$emission,
                                            current$converted_beta,
                                            current$mixing_coefficient,
                                            current$observations,
                                            params$n_group,
                                            params$n_states)
    log_likelihoods[, iter] <<- output[[1]]
    current$latent_sequences <<- output[[2]]
    current$z <<- output[[3]]
    samples$z[, iter] <<- current$z
    if (iter > params$warmup) {
      posterior$latent_sequences <<- posterior$latent_sequences + current$latent_sequences
    }
  }
  update_posterior <- function(select) {
    posterior$latent_sequences <<- posterior$latent_sequences / (params$n_iter - params$warmup)
    posterior$z <<- rowMeans(samples$z[, select])
  }
  if (training) {
    update_samples(iter)
  } else {
    update_posterior()
  }
}
