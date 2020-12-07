#' Fit Data to the Gompertz Distribution
#'
#' @description Fit mortality data to the Gompertz distribution using
#' nonlinear least squares. The algorithm used here is the Levenberg-Marquardt
#' algorithm. See \url{https://doi.org/10.1103/PhysRevE.83.036701} for more info.
#'
#' @param x vector of ages
#' @param y vector of number of deaths at age x
#' @param parstart vector of starting values for parameters \code{eta} and \code{b} respectively
#' @param eps stopping criterion
#' @param lambda damping parameter
#' @param lamup multiplied with lambda to decrease step size if the previous step is not accepted
#' @param lamdown lambda is divided by this to increase step size if the previous step is accepted
#'
#' @return List containing the vector of estimates for \code{eta} and \code{b} and the MSE at convergence.
#' @export
#'
#' @examples
#' x = seq(0, 5, by = 0.2)
#' y = abs(dgomp(x, eta = 0.1, b = 1) + rnorm(length(x), mean = 0, sd = 0.01))
#' gomp_fit(x, y)
gomp_fit <- function(x, y, parstart = NULL, eps = 0.0001, lambda = 1, lamup = 1.1, lamdown = 1.5){

  # Compatibility checks
  if(any(x < 0)){
    stop("x must be >= 0")
  }
  if(any(y < 0)){
    stop("y must be >= 0")
  }
  if(is.null(parstart)){
    parstart = c(1, 1)
  }else if(any(parstart) <= 0){
    stop("eta and b must be > 0")
  }

  # lamup and lamdown are used as part of our implementation of the Levenberg-Marquardt algorithm.
  # If a step is accepted, the step size is increased by multiplying lambda with lamup. If it is
  # not accepted, the step size is decreased by dividing lambda by lamdown, and another iteration
  # using the previously used parameter values (not the update) is used. The default values for
  # lamup and lamdown correspond to the "delayed gratification" method (Transtrum Machta Sethna 2011)
  # which speeds up the rate of convergence near the solution in exchange for slow progress toward the
  # beginning.

  if(lamup <= 1 | lamdown <= 1){
    stop("lamup/lamdown must be > 1")
  }

  out = gomp_fit_c(x, y, parstart = parstart, eps = eps, lambda = lambda, lamup = lamup, lamdown = lamdown)

  return(out)
}

