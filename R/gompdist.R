#' The Gompertz Distribution
#'
#' @description Density, distribution, and quantile functions along with random number
#' generation for the Gompertz distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param eta shape parameter
#' @param b scale parameter
#'
#' @return \code{dgomp} gives the density function, \code{pgomp} gives the distribution
#' function, \code{qgomp} gives the quantile function, and \code{rgomp} generates random values
#'
#' @examples
#' x = seq(0, 5, by = 0.1)
#' dgomp(x, eta = 0.2)
#' p <- 0:10 * 1 / 10
#' pgomp(qgomp(p, eta = 2), eta = 2)
#' rgomp(10, eta = 0.2)
#' @name gomp
NULL
#> NULL

#' @rdname gomp
#' @export
dgomp = function(x, eta, b = 1){
  if(eta <= 0){
    stop("Please enter a positive value for eta")
  }
  if(b <= 0){
    stop("Please enter a positive value for b")
  }
  if(any(x < 0)){
    stop("All x values must be >= 0")
  }
  return(b * eta * exp(eta + b * x - eta * exp(b * x)))
}

#' @rdname gomp
#' @export
pgomp = function(q, eta, b = 1){
  if(eta <= 0){
    stop("Please enter a positive value for eta")
  }
  if(b <= 0){
    stop("Please enter a positive value for b")
  }
  if(any(q < 0)){
    stop("All values of q must be >= 0")
  }
  return(1 - exp(-eta * (exp(b * q) - 1)))
}

#' @rdname gomp
#' @export
qgomp = function(p, eta, b = 1){
  if(eta <= 0){
    stop("Please enter a positive value for eta")
  }
  if(b <= 0){
    stop("Please enter a positive value for b")
  }
  if(any(p < 0) | any(p > 1)){
    stop("Probabilities must be between 0 and 1")
  }
  return((1 / b) * log(1 - (log(1 - p) / eta)))
}

#' @rdname gomp
#' @export
rgomp = function(n, eta, b = 1){
  if(eta <= 0){
    stop("Please enter a positive value for eta")
  }
  if(b <= 0){
    stop("Please enter a positive value for b")
  }
  u = runif(n)
  return((1 / b) * log(1 - (log(1 - u) / eta)))
}

