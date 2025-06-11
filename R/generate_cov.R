#' Generate an AR(1) correlation matrix
#'
#' Creates a square correlation matrix with an Autoregressive(1) structure.
#'
#' @param dim The dimension of the square matrix.
#' @param rho The AR(1) correlation parameter, a value in (-1, 1).
#'
#' @return A `dim`-by-`dim` matrix with the AR(1) structure.
#'
#' @export
#'
#' @importFrom stats toeplitz
generate_AR1 <- function(dim, rho) {
  #rho^abs(outer(1:p, 1:p, function(i, j) i - j))
  stats::toeplitz(rho^(0:(dim - 1)))
}

#' Generate a correlation matrix
#'
#' A helper function to generate a \eqn{K}-by-\eqn{K} correlation matrix, defaulting to an
#' AR(1) structure if a specific matrix `R` is not provided.
#'
#' @param rho The AR(1) correlation parameter, a value in (-1, 1). Used if `R` is `NULL`.
#' @param K The target dimension for the matrix. Required if `R` is `NULL`.
#' @param R An optional, pre-specified correlation matrix. If provided, `rho` and `K` are ignored.
#'
#' @return A \eqn{K}-by-\eqn{K} correlation matrix.
#'
#' @export
generate_R <- function(rho = NULL, K = NULL, R = NULL) {
  if (is.null(R)) {
    if (is.null(rho) || is.null(K)) {
      stop("`rho` and `K` must be specified when `R` is NULL.")
    }
    if (rho <= -1 || rho >= 1) {
      stop("AR(1) parameter `rho` must be in (-1, 1).")
    }
    if (!isTRUE(K > 1)) {
      stop("Dimension `K` must be an integer >= 2.")
    }
    R <- generate_AR1(dim = K, rho = rho)
  }
  return(R)
}

#' Generate a covariance matrix with a desired condition number ratio relative to a
#' correlation matrix
#'
#' Generates a \eqn{K}-by-\eqn{K} covariance matrix `Sigma` given a correlation
#' matrix `R` and a target condition number ratio `kappa.ratio` defines as the
#' ratio of spectral condition numbers \eqn{\kappa(\Sigma)/\kappa(R)}, where \eqn{\Sigma}
#' is the desired covariance matrix and \eqn{R} is the supplied correlation matrix.
#'
#' @param kappa.ratio The target condition number ratio, must be `kappa.ratio >= 1`.
#' @param R A \eqn{K}-by-\eqn{K} strictly positive definite correlation matrix.
#' @param q The number of "strong" (high relative variance) features. These features have
#'    their variance set to 1. The remaining `K - q` features have variance `s < 1`, where
#'    `s` is solved for numerically in order to satisfy the condition number ratio. Must
#'    be an integer `1 <= q < K`.
#' @param .tol.posdef Tolerance for checking whether the correlation matrix `R` is positive definite.
#' @param .tol.optimize Tolerance for the `optimize` function call.
#'
#' @return A list containing the generated `Sigma`, the input `R`, the condition
#'    numbers of both, target condition number ratio, and the final observed condition
#'    number ratio given the numerical optimization procedure.
#'
#' @export
#'
#' @importFrom stats optimize
#' @importFrom matrixcalc is.square.matrix is.positive.definite
generate_cov <- function(kappa.ratio, R, q = 1,
                         .tol.posdef = sqrt(.Machine$double.eps),
                         .tol.optimize = sqrt(.Machine$double.eps)) {
  K <- NROW(R)

  if (!is.numeric(q) || length(q) != 1 || q < 1 || q >= K) {
    stop("Number of strong features `q` must be an integer 1 <= q < NROW(R).")
  }
  q <- as.integer(q)
  if (!matrixcalc::is.square.matrix(R) || max(abs(diag(R) - 1)) > .Machine$double.eps) {
    stop("`R` must be a square correlation matrix with diagonal entries equal to 1.")
  }
  if (!matrixcalc::is.positive.definite(R, tol = .tol.posdef)) {
    stop("`R` must be positive definite.")
  }
  if (kappa.ratio < 1) {
    stop("Condition number ratio `kappa.ratio` must be >= 1.")
  }

  #--- internal helpers
  .make_SIG <- function(s, R.mat, q.val) {
    d.vec <- rep(c(1, s), times = c(q.val, NROW(R.mat) - q.val))
    R.mat * tcrossprod(d.vec)
  }
  .objective_SIG <- function(s, kappa.ratio.val, R.mat, q.val, kappa.R.val) {
    SIG <- .make_SIG(s = s, R.mat = R.mat, q.val = q.val)
    abs(kappa(SIG, exact = TRUE) - kappa.ratio.val * kappa.R.val)
  }
  #--- end helpers

  kappa.R <- kappa(R, exact = TRUE)

  opt <- stats::optimize(
    .objective_SIG,
    lower = .tol.optimize,
    upper = 1,
    kappa.ratio.val = kappa.ratio,
    R.mat = R,
    q.val = q,
    kappa.R.val = kappa.R,
    tol = .tol.optimize
  )

  Sigma <- .make_SIG(s = opt$minimum, R.mat = R, q.val = q)
  if (!isSymmetric(Sigma, tol = .tol.posdef)) Sigma <- (Sigma + t(Sigma))/2
  kappa.Sigma <- kappa(Sigma, exact = TRUE)

  list(
    R = R,
    Sigma = Sigma,
    kappa.R = kappa.R,
    kappa.Sigma = kappa.Sigma,
    kappa.ratio.tar = kappa.ratio,
    kappa.ratio.obs = kappa.Sigma / kappa.R
  )
}
