## armor
armor assists those interested in studying age-related mortality by
providing tools for generating data from and estimating the parameters of the probability
distributions commonly used in aging research.

# Installation
    # Install directly from GitHub:
    # install.packages("devtools")
    devtools::install_github("donrturner/armor")

# Introduction
Barring premature death, senescence (i.e., biological aging) is a condition which one hundred percent of people worldwide are susceptible to. It is well known that the risk of death increases as one ages. Though we have in recent years obtained a far greater understanding of aging, our understanding is far from complete. This package aims to provide tools to those interested in studying senescence, specifically concerning the relationship between age and mortality.

Some tools used to model age-related mortality are the Gompertz and Gompertz-Makeham distributions (Wilson, 1994). This package provides tools for simulating observations from these distributions, as well as for parameter estimation.

# Simulation
Functions are provided for simulating observations from the Gompertz and Gompertz-Makeham distributions. Included are functions for obtaining values from the density, distribution, and quantile functions, as well as functions for generating random observations.

## Examples
```
# Gompertz distribution
x = seq(0, 5, by = 0.1)
dgomp(x, eta = 0.2)
p <- 0:10 * 1 / 10
pgomp(qgomp(p, eta = 2), eta = 2)
rgomp(10, eta = 0.2)

# Gompertz-Makeham Distribution
x = seq(0.5, 5, by = 0.1)
dgompm(x, 1, 1, 1)
p <- seq(0, .99, by = 0.1)
pgompm(qgompm(p, .5, .6, .7), .5, .6, .7)
rgompm(10, 1, 2, 3)
```

# Estimation
A function for fitting data to the Gompertz density function is provided. Data is fit via the Levenberg-Marquardt algorithm with delayed gratification for nonlinear least squares. Click [here](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.83.036701#fulltext) for more info.

## Example
```
x = seq(0, 5, by = 0.2)
y = abs(dgomp(x, eta = 0.1, b = 1) + rnorm(length(x), mean = 0, sd = 0.01))
gomp_fit(x, y)
```

# References
Wilson, D. L. (1994). The analysis of survival (mortality) data: Fitting Gompertz, Weibull, and logistic functions. *Mechanisms of Ageing and Development*, 15-33.

