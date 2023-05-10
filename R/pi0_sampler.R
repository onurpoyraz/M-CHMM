pi0_sampler <- function(iter, training = TRUE, select = NULL) {
  update_samples <- function(iter) {
    for (group in 1:params$n_group) {
      if (iter == 1){
        samples$pi0[[group]] <<- list()
        current$pi0[[group]] <<- list()
      }
      for (chain in 1:params$n_chains) {
        if (iter == 1) {
          samples$pi0[[group]][[chain]] <<- array(as.numeric(NA), c(params$n_states, params$n_iter))
          current$pi0[[group]][[chain]] <<- prior$pi0[[chain]] / sum(prior$pi0[[chain]])
        } else {
          current$pi0[[group]][[chain]] <<- sample_pi0(prior$pi0[[chain]],
                                                       current$latent_sequences,
                                                       current$z,
                                                       params$n_states,
                                                       as.integer(group - 1),
                                                       as.integer(chain - 1))
        }
        samples$pi0[[group]][[chain]][ , iter] <<- current$pi0[[group]][[chain]]
      }
    }
  }
  update_posterior <- function(select) {
    for (group in 1:params$n_group) {
      posterior$pi0[[group]] <<- list()
      for (chain in 1:params$n_chains) {
        posterior$pi0[[group]][[chain]] <<- rowMeans(samples$pi0[[group]][[chain]][ , select])
      }
      names(samples$pi0[[group]]) <<- params$chain_names
      names(posterior$pi0[[group]]) <<- params$chain_names
    }
  }
  if (training) {
    update_samples(iter)
  } else {
    update_posterior(select)
  }
}
