#' Solve for OLS coefficients on a subset of data
#'
#' A stable OLS solver using QR decomposition. This is an internal helper
#' function for computing coefficients within a child node.
#'
#' @param indices A logical or integer vector of indices specifying the subset of
#'    data (the child node) to use.
#' @param w.mat The full matrix of regressors.
#' @param y.vec The full vector of outcomes.
#' @param intercept Logical, whether to include an intercept in the regression.
#'    Defaults to `intercept = TRUE`.
#'
#' @return A numeric vector of the solved OLS coefficients.
#'
#' @keywords internal
solve_OLS <- function(indices, w.mat, y.vec, intercept = TRUE) {
  w.mat <- w.mat[indices, , drop = FALSE]
  y.vec <- y.vec[indices]
  if (intercept) w.mat <- cbind(1, w.mat)
  drop(qr.solve(w.mat, y.vec))
}

#' Compute empirical solution for the local parameter estimates under a varying-coefficient model
#' fit across a candidate split
#'
#' A wrapper function that computes child solutions \eqn{\hat\theta_{C_1}} and \eqn{\hat\theta_{C_2}}
#' for the optimal local parameters fit over the data in \eqn{C_1} and \eqn{C_2}, respectively.
#'
#' @param C1 A vector of indices (logical or integer) for the first child node. See [solve_OLS()].
#' @param C2 A vector of indices (logical or integer) for the second child node. See [solve_OLS()].
#' @inheritParams solve_OLS
#' @param labels An optional character vector to name the resulting coefficients.
#' @param ... Additional arguments passed to [solve_OLS()].
#'
#' @return A list containing the named coefficient vectors for C1 and C2.
#'
#' @keywords internal
compute_theta_hat <- function(C1, C2, w.mat, y.vec, labels = NULL, ...) {
  lapply(list("C1" = C1, "C2" = C2), function(indices) {
    # The first coefficient from solve_OLS is the intercept, which we drop to
    # keep only the theta components for the criteria calculations.
    res <- solve_OLS(indices = indices, w.mat = w.mat, y.vec = y.vec, ...)[-1]
    if (!is.null(labels)) names(res) <- labels
    return(res)
  })
}
