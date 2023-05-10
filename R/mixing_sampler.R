mixing_sampler <- function(iter, training = TRUE, select = NULL) {
  update_samples <- function(iter) {
    if (iter == 1) {
      samples$mixing_coefficient <<- array(as.numeric(NA), c(params$n_group, params$n_iter))
      current$mixing_coefficient <<- prior$mixing_coefficient / sum(prior$mixing_coefficient)
    } else {
      current$mixing_coefficient <<- sample_mixing_coefficient(prior$mixing_coefficient,
                                                               current$z)
    }
    samples$mixing_coefficient[, iter] <<- current$mixing_coefficient
  }
  update_posterior <- function(select) {
    posterior$mixing_coefficient <<- rowMeans(samples$mixing_coefficient[, select])
  }
  if (training) {
    update_samples(iter)
  } else {
    update_posterior(select)
  }
}
