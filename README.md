## armor
armor assists those interested in studying age-related mortality by
providing tools for generating data from and fitting models to the probability
distributions commonly used in aging research.

# Installation
    # Install directly from GitHub:
    # install.packages("devtools")
    devtools::install_github("donrturner/armor")

# Introduction
Barring premature death, senescence (i.e., biological aging) is a condition which one hundred percent of people worldwide are susceptible to. It is well known that the risk of death increases as one ages. Though we have in recent years obtained a far greater understanding of aging, our understanding is far from complete. This package aims to provide tools to those interested in studying senescence, specifically concerning the relationship between age and mortality.

Some tools used to model age-related mortality are the Gompertz distribution and Gompertz-Makeham Law of Mortality (Wilson, 1994). Another common distribution used in longevity analysis is the Logistic–Makeham distribution (Pletcher, 2000). This package provides tools for simulating observations from these distributions, as well as for creating fitted models.

# References
Pletcher, S. D. (2000). Why Do Life Spans Differ? Partitioning Mean Longevity Differences in Terms of Age-Specific Mortality Parameters. *The Journals of Gerontology: Series A*, B381–B389.

Wilson, D. L. (1994). The analysis of survival (mortality) data: Fitting Gompertz, Weibull, and logistic functions. *Mechanisms of Ageing and Development*, 15-33.

