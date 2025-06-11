# criteval: A basic simulation pipeline for evaluating GRF-style splitting criteria

This package provides a simple, encapsulated framework for running some basic simulations that evaluate different splitting criteria used in tree-based recursive partitioning, particularly those in the style of Generalized Random Forests (GRF).

## Installation

Should anyone else be interested in playing with this, you can install the development version of `criteval` from GitHub with:
```R
# install.packages("remotes")
remotes::install_github("dfleis/criteval")
```

## TODO

* Multivariate case ($K > 2$, $p > 1$).

* (Readme) Add our GRF-FPT reference. Provide a more complete summary and review of the problem for better context.

## Scope and intention

* **Data generation:** Includes tools to generate synthetic data following varying-coefficient models where the correlation and conditioning of the regressors are precisely controlled.
The underlying varying-coefficient model looks like:
$$
\mathbb E[Y \mid X = x] = \nu^*(x) + W^\top \theta^*(x), 
$$
where $W = (W_1,\ldots, W_K)^\top \in \mathbb R^K$ denotes a set of primary regressors, $Y \in \mathbb R$ denotes a scalar outcome, and $X \in \mathcal X = [0,1]^p$ denotes auxiliary covariates such that the model can be described as *conditionally linear* given $X$.
Here, the $\theta^*(x) = (\theta_1^*(x),\ldots,\theta_K^*(x))$ denote the target effect functions with each component function $\theta_k^*(x)$ denoting the effect of the corresponding regressor $W_k$ on the outcome $Y$ local to the covariates at $X = x$.

* **Criteria evaluation:** The main simulation loop is designed to scan across all potential binary splits in the initial dataset and calculates the value of multiple criteria at each candidate split.

## References

Susan Athey, Julie Tibshirani and Stefan Wager.
<b>Generalized Random Forests.</b> <i>Annals of Statistics</i>, 47(2), 2019.
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>,
<a href="https://arxiv.org/abs/1610.01271">arxiv</a>]
