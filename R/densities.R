dlaplace <- function(data, mu, scale, log=TRUE){
  x <- - abs(data - mu) / scale
  y <- log(0.5) - log(scale)
  log_density <- x + y
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}
dhorseshoe <- function(data, mu, scale, log=TRUE) {
  x <- ((data - mu) / scale)^2 / 2
  g <- 0.5614594835668851  # exp(-0.5772156649015328606)
  b <- 1.0420764938351215  # sqrt(2 * (1-g) / (g * (2-g)))
  h_inf <- 1.0801359952503342  #  (1-g)*(g*g-6*g+12) / (3*g * (2-g)**2 * b)
  q <- (20. / 47.) * x^1.0919284281983377
  h <- 1. / (1 + x^(1.5)) + (h_inf * q) / (1 + q)
  f <- 0.5 * log(2 * pi^3) + log(g * scale)
  log_density <- log(log(1 + g / (x + 1e-12) - (1 - g) / (h + b * x)^2)) - log(1 + (1 - g) / g * exp(-x / (1 - g))) - f
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}
prior_density <- function(beta_values, distribution, sd) {
  if (length(beta_values) == 0) {
    return (0)
  } else if (distribution == "Gaussian") {
    return (sum(dnorm(beta_values, 0, sd, log = TRUE)))
  } else if (distribution == "Laplace") {
    return (sum(dlaplace(beta_values, 0, sd, log = TRUE)))
  } else if (distribution == "Horseshoe") {
    return (sum(dhorseshoe(beta_values, 0, sd, log = TRUE)))
  }
}