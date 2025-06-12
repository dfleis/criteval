# criteval: A toy simulation pipeline for evaluating GRF-style splitting criteria for identifying parameter heterogeneity under varying-coefficient models

This package provides a simple, encapsulated framework for running some basic simulations that evaluate different splitting criteria used in tree-based recursive partitioning, particularly those in the style of Generalized Random Forests (GRF).
This framework looks at the behaviour of several splitting criteria designed to identify heterogeneity in the data with respect to estimates of some underlying effect functions $\theta^*(x)$ in under a varying-coefficient model (VCM):

$$
\mathbb E[Y \mid X = x] = \nu^\star (x) + W^\top \theta^\star (x), 
$$

Here, $W = (W_1,\ldots, W_K)^\top \in \mathbb R^K$ denotes a set of primary regressors, $Y \in \mathbb R$ denotes a scalar outcome, and $X \in \mathcal X = [0,1]^p$ denotes auxiliary covariates such that the model can be described as *conditionally linear* given $X$.
Here, the $\theta^\star(x) = (\theta_1^\star(x),\ldots,\theta_K^\star(x))$ denote the target effect functions with each component function $\theta_k^*(x)$ denoting the effect of the corresponding regressor $W_k$ on the outcome $Y$ local to the covariates at $X = x$.

## Installation

Should anyone else be interested in playing with this, you can install the development version of `criteval` from GitHub with:
```R
# install.packages("remotes")
remotes::install_github("dfleis/criteval")
```

## TODO

* Implement the special case for bivariate $W$ which has a closed-form solution for the covariance/condition number.

* Change bivariate case to be written in terms of the condition number ratio $r := \kappa(\Sigma)/\kappa(R) \geq 1$. 
For two random variables with correlation coefficient $\rho$, one has $\kappa(R) = (1 + |\rho|)/(1 - |\rho|)$. This is equivalent to
the current method of directly supplying a condition number for $\kappa(\Sigma)$, but is more straightforward to work with since all one has to worry
about is whether the ratio $r \geq 1$.

* More descriptive variable names, particularly in `create_data_generator.R`

* (Readme) Add our GRF-FPT reference. Provide a more complete summary and review of the problem for better context.

## Scope and intention

* **Data generation:** Includes tools to generate synthetic data following varying-coefficient models where the correlation and conditioning of the regressors are precisely controlled.

* **Criteria evaluation:** The main simulation loop is designed to scan across all potential binary splits in the initial dataset and calculates the value of multiple criteria at each candidate split.

## References

Susan Athey, Julie Tibshirani and Stefan Wager.
<b>Generalized Random Forests.</b> <i>Annals of Statistics</i>, 47(2), 2019.
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>,
<a href="https://arxiv.org/abs/1610.01271">arxiv</a>]
