#' Create a data generator factory
#'
#' This function acts as a factory method to produce a data generator object.
#' It takes all the parameters that define a data-generating process (DGP) and
#' returns an object that can be used to repeatedly generate datasets from that
#' DGP.
#'
#' The regressors `W` are sampled from a multivariate Gaussian whose covariance
#' structure is determined by the provided correlation structure (either through
#' the AR(1) correlation parameter `rho` or via an explicit correlation matrix `R`)
#' as well as target condition number ratio `kappa.ratio`. The condition number
#' ratio is the ratio of the condition numbers for the desired covariance matrix
#' and generated/supplied correlation matrix.
#'
#' The auxiliary covariates `X` are generated from a Gaussian copula with an
#' AR(1) correlation structure defined by `rho.X`. When `rho.X = 0` (the default),
#' the `X` covariates are independent Uniform(0,1).
#'
#' This factory ensures reproducibility. Providing the same `seed` will yield
#' the same random components, as the main seed is used to derive a stable
#' set of component-specific seeds.
#'
#' @inheritParams generate_cov
#' @param rho.W The AR(1) correlation parameter for the `W` regressors. See `rho` in [generate_corr()].
#'    Ignored if a correlation matrix `R.W` is provided.
#' @param R.W An optional, pre-specified correlation matrix for `W`. See `R` in [generate_cov()].
#' @param q.W The number of "strong" (high relative variance) features in `W`. See `q` in [generate_cov()].
#' @param mu.W An optional mean vector for the regressors `W`. Defaults to a zero vector.
#' @param p.X The dimension of the auxiliary covariates `X`.
#' @param rho.X The AR(1) correlation for the latent Gaussian variables used to
#'    generate the `X` covariates via a copula. Defaults to `rho.X = 0`, which
#'    corresponds to independent Uniform(0,1) covariates.
#' @param sigma.eps The standard deviation of the Gaussian noise term epsilon.
#' @param nu_FUN A function that takes an `n`-by-`p.X` matrix `X` and returns the
#'    `n`-dimensional nuisance vector `nu_FUN(X)`.
#' @param theta_FUN_list A list of \eqn{K} functions, where each function takes the
#'    `n`-by-`p.X` matrix `X` and returns an `n`-dimensional vector for one component
#'    of `theta_FUN_list[[k]](X)`.
#' @param ... Additional arguments passed to [generate_cov()].
#'
#' @return An object of class `criteval_generator`. This is a list containing
#'    generator functions (`$generate`, `$X`, `$W`) and key parameters (`$Sigma`,
#'    `$theta_FUN`, etc.) that define the DGP.
#'
#' @examples
#' \dontrun{
#' # Define the DGP parameters and create the generator
#' theta_funs <- list(
#'   theta1 = function(x) plogis((x[,1] - 0.45), scale = 1/10),
#'   theta2 = function(x) plogis(-(x[,1] - 0.55), scale = 1/10)
#' )
#' nu_fun <- function(x) { 0 * x[,1] }
#'
#' gen <- data_generator(
#'   kappa.ratio = 1.5, # W: Target condition number ratio kappa(Sigma.W)/kappa(R.W)
#'   rho.W = 0.6,       # W: AR(1) correlation parameter for W
#'   q.W = 1,           # W: Number of relative high variance regressors W
#'   p.X = 2,           # X: Auxiliary covariate dimension p
#'   rho.X = 0.8,       # X: AR(1) correlation parameter for X
#'   sigma.eps = 0.1,   # Gaussian noise standard deviation
#'   nu_FUN = nu_fun,   # VCM intercept (nuisance function)
#'   theta_FUN_list = theta_funs # Underlying target effect functions
#' )
#'
#' # Use the generator to create datasets
#' data <- gen$generate(n = 5, seed = 123)
#'
#' # # Reproducible data
#' # data$X
#' # gen$X(n = 5, seed = 123)
#'
#' # data$W
#' # gen$W(n = 5, seed = 123)
#' }
#'
#' @export
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats pnorm rnorm
data_generator <- function(
    kappa.ratio = 1, rho.W = NULL, R.W = NULL, q.W = 1, mu.W = NULL,
    p.X = 1, rho.X = 0, sigma.eps = 1, theta_FUN_list, nu_FUN, ...) {

  #--- One-time setup based on DGP parameters
  .make_theta_FUN <- function(FUN_list) {
    if (!is.list(FUN_list) || !all(sapply(FUN_list, is.function))) {
      stop("`FUN_list` must be a list of functions.")
    }
    function(x) sapply(FUN_list, function(f) f(x))
  }
  theta_FUN <- .make_theta_FUN(theta_FUN_list)
  K.W <- length(theta_FUN_list)

  # Setup for W regressors
  R.W <- generate_corr(R = R.W, rho = rho.W, K = K.W)
  cov.obj.W <- generate_cov(kappa.ratio = kappa.ratio, R = R.W, q = q.W, ...)
  Sigma.W <- cov.obj.W$Sigma
  if (is.null(mu.W)) mu.W <- rep(0, K.W)
  if (length(mu.W) != nrow(Sigma.W)) {
    stop(
      "Mean vector `mu.W` and covariance matrix `Sigma.W` have non-conforming size:",
      sprintf(
        "\n\tlength(mu.W) = %i, NROW(Sigma.W) = %i, length(theta_FUN_list) = %i",
        length(mu.W), NROW(R.W), length(theta_FUN_list)
      )
    )
  }

  # Setup for X covariates (latent correlation)
  Sigma.X.latent <- generate_AR1(dim = p.X, rho = rho.X)

  #--- Define component generator functions (as closures)
  generate_X_component <- function(n, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    Z.latent <- mvtnorm::rmvnorm(n, mean = rep(0, p.X), sigma = Sigma.X.latent)
    X <- apply(Z.latent, 2, stats::pnorm)
    colnames(X) <- paste0("X", 1:p.X)
    structure(X, Sigma.latent = Sigma.X.latent)
  }

  generate_W_component <- function(n, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    W <- mvtnorm::rmvnorm(n = n, mean = mu.W, sigma = Sigma.W)
    colnames(W) <- paste0("W", 1:K.W)
    structure(W, mu = mu.W, R = R.W, Sigma = Sigma.W)
  }

  generate_full_dataset <- function(n, seed = NULL) {
    if (!is.null(seed)) {
      set.seed(seed)
      sub.seeds <- sample.int(.Machine$integer.max, 3)
      seed.X <- sub.seeds[1]
      seed.W <- sub.seeds[2]
      seed.eps <- sub.seeds[3]
    } else {
      seed.X <- seed.W <- seed.eps <- NULL
    }

    X <- generate_X_component(n, seed = seed.X)
    W <- generate_W_component(n, seed = seed.W)

    if (!is.null(seed.eps)) set.seed(seed.eps)
    eps <- stats::rnorm(n, 0, sigma.eps)

    nuX <- nu_FUN(X)
    thetaX <- theta_FUN(X)
    Y <- nuX + rowSums(thetaX * W) + eps
    list(X = X, Y = Y, W = W, nuX = nuX, thetaX = thetaX)
  }

  #--- Assemble and return the generator object
  generator <- list(
    "generate" = generate_full_dataset,
    "X" = function(n, seed = NULL) {
      seed_X <- if (!is.null(seed)) { set.seed(seed); sample.int(.Machine$integer.max, 1) } else { NULL }
      generate_X_component(n, seed = seed_X)
    },
    "W" = function(n, seed = NULL) {
      seed_W <- if (!is.null(seed)) { set.seed(seed); sample.int(.Machine$integer.max, 3)[2] } else { NULL }
      generate_W_component(n, seed = seed_W)
    },
    "theta_FUN" = theta_FUN,
    "nu_FUN" = nu_FUN,
    "Sigma.W" = Sigma.W,
    "R.W" = R.W,
    "Sigma.X.latent" = Sigma.X.latent,
    "params" = list(
      K.W         = K.W,
      kappa.ratio = kappa.ratio,
      rho.W       = rho.W,
      q.W         = q.W,
      mu.W        = mu.W,
      p.X         = p.X,
      rho.X       = rho.X,
      sigma.eps   = sigma.eps
    )
  )

  class(generator) <- "criteval_data_generator"
  return(generator)
}
