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

