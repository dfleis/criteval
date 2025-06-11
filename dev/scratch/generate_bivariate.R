#' Generate bivariate Gaussian samples
#'
#' Generates a 2-column matrix containing samples of a bivariate Gaussian vector
#' with a specified underlying correlation coefficient and a specified spectral
#' condition for its covariance matrix.
#'
#' @param n The number of observations to generate.
#' @param corr The desired correlation between the two variables, a value in (-1, 1).
#' @param kappa.cov The desired spectral condition number of the covariance matrix. If `NULL`,
#'    it is set to the minimum possible value for the given `corr`. See [calc_kappa_lb_2d()].
#' @param mu A numeric vector of length 2 for the mean of the variables. Defaults to `mu = c(0, 0)`.
#' @param .tol A small tolerance for numerical comparisons.
#'
#' @return An n-by-2 matrix of regressors. The matrix has attributes
#'   `mu`, `Sigma`, `corr`, and `kappa.cov` attached.
#'
#' @export
generate_bivariate <- function(n, corr, kappa.cov = NULL, mu = c(0, 0), .tol = sqrt(.Machine$double.eps)) {
  if (any(!is.numeric(n), n <= 0)) stop("Sample size `n` must be a positive integer.")
  if (any(!is.numeric(mu), length(mu) != 2)) stop("Population mean `mu` must be a numeric vector of length 2.")
  kappa.info <- validate_corr_kappa_2d(corr = corr, kappa.cov = kappa.cov, .tol = .tol)

  #--- Calculate variance ratio
  if (abs(kappa.info$kappa.cov - kappa.info$kappa.lb) < .tol) {
    r <- 1 # If very close to lower bound, set variance ratio r = 1
  } else {
    # Solve quadratic for variance ratio
    A <- kappa.info$kappa.cov
    B <- corr^2 * (kappa.info$kappa.cov + 1)^2 - (kappa.info$kappa.cov^2 + 1)
    d <- B^2 - 4 * A^2

    if (all(d < 0, abs(d) < .tol)) {
      warning("Numerical issue: Negative discriminant. Discriminant to zero since it was within tolerance `.tol`. Check inputs.")
      d <- 0
    } else if (d < 0) {
      # Can this happen under a feasible condition number?
      stop("Negative discriminant that exceeds the tolerance.")
    }
    r <- (-B + sqrt(d)) / (2 * A)
  }

  #--- Construct covariance matrix
  sigma1 <- 1
  sigma2 <- 1/sqrt(r)
  Sigma <- matrix(
    c(
      sigma1^2, kappa.info$corr * sigma1 * sigma2,
      kappa.info$corr * sigma1 * sigma2, sigma2^2
    ),
    nrow = 2,
    byrow = TRUE
  )

  W <- mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)
  structure(W, mu = mu, Sigma = Sigma, corr = kappa.info$corr, kappa.cov = kappa.info$kappa.cov)
}


#' Validate the feasibility between a correlation coefficient and a condition number
#'
#' Checks whether the desired spectral condition number of some 2x2 covariance matrix is
#' feasible given the supplied correlation coefficient.
#'
#' @param corr The correlation coefficient between two random variables.
#' @param kappa.cov The spectral condition number of a 2x2 covariance matrix for random variables
#'    with correlation `corr`.
#' @param .tol A small tolerance for numerical comparisons.
#'
#' @return Returns a named list containing the original input correlation `corr`,
#'    the validated condition number `kappa.cov`, and the lower bound `kappa.lb` for
#'    the condition number of the covariance between two variables with correlation `corr`.
#'    See [calc_kappa_lb_2d()].
#'
#' @export
validate_corr_kappa_2d <- function(corr, kappa.cov = NULL, .tol = sqrt(.Machine$double.eps)) {
  if (any(!is.numeric(corr), abs(corr) >= 1)) stop("Correlation `corr` must satisfy -1 < corr < 1.")
  if (is.null(kappa.cov)) kappa.cov <- calc_kappa_lb_2d(corr)
  if (any(!is.numeric(kappa.cov), kappa.cov < 1)) stop("Covariance condition number `kappa.cov` must satisfy kappa.cov >= 1.")
  kappa.lb <- calc_kappa_lb_2d(corr)

  if (abs(kappa.cov - kappa.lb) < .tol) kappa.cov <- kappa.lb

  if (kappa.cov < kappa.lb) {
    stop(
      "Condition number `kappa.cov` must satisfy must satisfy kappa.cov >= (1 + |corr|)/(1 - |corr|).\n",
      "This lower bound represents the smallest possible condition number of the covariance given two variables with correlation `corr`.\n",
      sprintf("\nThe target correlation `corr = %s` implies the lower bound:\n", corr),
      sprintf("\t(1 + |corr|)/(1 - |corr|) = %0.4f\n", round(kappa.lb, 4)),
      sprintf("but supplied condition number is `kappa.cov = %s` with an gap of:\n", kappa.cov),
      sprintf("\tkappa.cov - kappa.lb = %0.2E", kappa.cov - kappa.lb)
    )
  }
  list(
    corr = corr,
    kappa.cov = kappa.cov,
    kappa.lb = kappa.lb
  )
}

#' Calculate the lower bound of the condition number for the covariance matrix
#' of two random variables with a given correlation
#'
#' The condition number of a covariance matrix \eqn{\Sigma = \text{Cov}(W)} for
#' a bivariate random vector \eqn{W = (W_1, W_2)} is bounded below by
#' \deqn{\kappa(\Sigma) \geq \frac{1 + |\rho|}{1 - |\rho|}}
#' where \eqn{\rho = \text{Corr}(W_1,W_2) \in (-1, 1)} is the underlying correlation
#' between the two component variables.
#'
#' @param corr The correlation coefficient, a numeric value between -1 and 1.
#'
#' @return The numeric lower bound of the spectral condition number.
#'
#' @export
calc_kappa_lb_2d <- function(corr) {
  (1 + abs(corr)) / (1 - abs(corr))
}
