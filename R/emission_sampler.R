emission_sampler <- function(iter, training = TRUE, select = NULL) {
  update_samples <- function(iter) {
    for (group in 1:params$n_group) {
      if (iter == 1){
        samples$emission[[group]] <<- list()
        current$emission[[group]] <<- list()
      }
      for (chain in 1:params$n_chains) {
        if (iter == 1) {
          samples$emission[[group]][[chain]] <<- array(as.numeric(NA), c(params$n_states, params$n_observations, params$n_iter))
          current$emission[[group]][[chain]] <<- prior$emission[[chain]] / rowSums(prior$emission[[chain]])
        } else if (iter > params$emission_lag) {
          current$emission[[group]][[chain]] <<- sample_emission(prior$emission[[chain]],
                                                                 current$latent_sequences,
                                                                 current$observations,
                                                                 current$z,
                                                                 params$n_states,
                                                                 params$n_observations, 
                                                                 as.integer(group - 1),
                                                                 as.integer(chain - 1))
        }
        samples$emission[[group]][[chain]][ , , iter] <<- current$emission[[group]][[chain]]
      }
    }
  }
  update_posterior <- function(select) {
    for (group in 1:params$n_group) {
      posterior$emission[[group]] <<- list()
      for (chain in 1:params$n_chains) {
        posterior$emission[[group]][[chain]] <<- apply(samples$emission[[group]][[chain]][ , , select], c(1,2), FUN=mean)
      }
      names(samples$emission[[group]]) <<- params$chain_names
      names(posterior$emission[[group]]) <<- params$chain_names
    }
  }
  if (training) {
    update_samples(iter)
  } else {
    update_posterior(select)
  }
}
