
#' The Gompertz-Makeham Distribution
#'
#' @description Density, distribution, and quantile functions along with random number
#' generation for the Gompertz-Makeham distribution.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param alpha,beta,lambda parameters
#'
#' @return \code{dgompm} gives the density function, \code{pgompm} gives the distribution
#' function, \code{qgompm} gives the quantile function, and \code{rgompm} generates random values
#'
#' @examples
#' x = seq(0.5, 5, by = 0.1)
#' dgompm(x, 1, 1, 1)
#' p <- 0:10 * 1 / 10
#' pgompm(qgompm(p, .5, .6, .7), .5, .6, .7)
#' rgompm(10, 1, 2, 3)
#' @name gompm
NULL
#> NULL

# Gompertz-Makeham density
#' @rdname gompm
#' @export
dgompm = function(x, alpha, beta, lambda){
  if(alpha <= 0){
    stop("Please enter a positive value for alpha")
  }
  if(beta <= 0){
    stop("Please enter a positive value for beta")
  }
  if(any(lambda <= 0)){
    stop("Please enter a positive value for lambda")
  }
  if(any(x <= 0)){
    stop("All x values must be > 0")
  }
  return((alpha * exp(beta * x) + lambda) * exp(-lambda * x - (alpha / beta) * (exp(beta * x) - 1)))
}

# Gompertz-Makeham distribution function
#' @rdname gompm
#' @export
pgompm = function(q, alpha, beta, lambda){
  if(alpha <= 0){
    stop("Please enter a positive value for alpha")
  }
  if(beta <= 0){
    stop("Please enter a positive value for beta")
  }
  if(any(lambda <= 0)){
    stop("Please enter a positive value for lambda")
  }
  if(any(q <= 0)){
    stop("All q values must be > 0")
  }
  return(1 - exp(-lambda * q - (alpha / beta) * (exp(beta * q) - 1)))
}

# Gompertz-Makeham quantile function
#' @rdname gompm
#' @export
qgompm = function(p, alpha, beta, lambda){
  if(alpha <= 0){
    stop("Please enter a positive value for alpha")
  }
  if(beta <= 0){
    stop("Please enter a positive value for beta")
  }
  if(any(lambda <= 0)){
    stop("Please enter a positive value for lambda")
  }
  if(any(p < 0) | any(p > 1)){
    stop("Probabilities must be between 0 and 1")
  }
  # The Lambert W function is required in the formula of the quantile function for a closed-form expression
  return(alpha / (beta * lambda) - (1 / lambda) * log(1 - p) - (1 / beta) * lambertW0(alpha * exp(alpha / lambda) * (1 - p)^(-beta / lambda) / lambda))
}

# Generate random values from the Gompertz-Makeham distribution
#' @rdname gompm
#' @export
rgompm = function(n, alpha, beta, lambda){
  if(alpha <= 0){
    stop("Please enter a positive value for alpha")
  }
  if(beta <= 0){
    stop("Please enter a positive value for beta")
  }
  if(any(lambda <= 0)){
    stop("Please enter a positive value for lambda")
  }
  u = runif(n)
  # The Lambert W function is required in the formula of the quantile function for a closed-form expression
  return(alpha / (beta * lambda) - (1 / lambda) * log(1 - u) - (1 / beta) * lambertW0(alpha * exp(alpha / lambda) * (1 - u)^(-beta / lambda) / lambda))
}


