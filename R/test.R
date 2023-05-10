fixed_pi0 <- function(iter, training = TRUE, select = NULL) {
  update_samples <- function(iter) {
    for (group in 1:params$n_group) {
      if (iter == 1) {
        current$pi0[[group]] <<- list()
        sample <- test$pi0[[group]]
        dim(sample) <- params$pi0_dim
        for (chain in 1:params$n_chains) {
          if (params$n_covariates == 0) {
            pi0 <- apply(sample, chain, FUN=sum)
            current$pi0[[group]][[chain]] <<- matrix(pi0 / sum(pi0), ncol=1)
          } else {
            pi0 <- apply(sample, chain, FUN=c)
            current$pi0[[group]][[chain]] <<- t(pi0 / rowSums(pi0))
          }
        }
      }
    }
  }
  update_posterior <- function(select) {
    for (group in 1:params$n_group) {
      posterior$pi0[[group]] <<- test$pi0[[group]]
    }
  }
  if (training) {
    update_samples(iter)
  } else {
    update_posterior(select)
  }
}


fixed_emission <- function(iter, training = TRUE, select = NULL) {
  update_samples <- function(iter) {
    for (group in 1:params$n_group) {
      if (iter == 1) {
        current$emission[[group]] <<- test$emission[[group]]
      }
    }
  }
  update_posterior <- function(select) {
    for (group in 1:params$n_group) {
      posterior$emission[[group]] <<- test$emission[[group]]
    }
  }
  if (training) {
    update_samples(iter)
  } else {
    update_posterior(select)
  }
}


fixed_beta <- function(iter, training = TRUE, select = NULL) {
  update_samples <- function(iter) {
    for (group in 1:params$n_group) {
      if (iter == 1) {
        current$beta[[group]] <<- test$beta[[group]]
        current$transition[[group]] <<- test$transition[[group]]
      }
    }
  }
  update_posterior <- function(select) {
    for (group in 1:params$n_group) {
      posterior$beta[[group]] <<- test$beta[[group]]
      posterior$transition[[group]] <<- test$transition[[group]]
    }
  }
  if (training) {
    update_samples(iter)
  } else {
    update_posterior(select)
  }
}

fixed_latent_sequence <- function (iter, training = TRUE) {
  update_samples <- function(iter) {
    if (iter == 1) {
      latent_sequences <<- test$latent_sequences
    }
  }
  update_posterior <- function(select) {
    posterior$latent_sequences <<- test$latent_sequences
  }
  if (training) {
    update_samples(iter)
  } else {
    update_posterior(select)
  }
}