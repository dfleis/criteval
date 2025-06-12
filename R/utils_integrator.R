#' Pre-computes inner integrals to create a lookup table
#'
#' This function performs expensive, high-dimensional integration upfront. For a
#' given candidate splitting axis/dimension \eqn{\ell \in \{1,\ldots,p\}} over the
#' covariate space, it integrates the target function \eqn{\theta^*(x)} over the
#' other \eqn{p - 1} dimensions at various discrete points along the \eqn{\ell}-th axis.
#' The result is a lookup table that can be queried quickly and used for interpolation
#' with faster, one-dimensional integration.
#'
#' Let's say we want to calculate the integral of \eqn{\theta^*(x)} over a hyperrectangle
#' (child node) \eqn{C = [0,t] \times [0,1]^{p-1}}, where we're assuming for illustration
#' that the candidate split is along the first dimension \eqn{\ell = 1}. The full integral is
#' \deqn{\int_{C} \theta^*(x)\,dx = \int_0^t \left( \int_{[0,1]^{p-1}}\theta^*(x_1,\ldots,x_p)\,dx_2\cdots\,dx_p \right)\,dx_1}
#' This breaks the problem into two parts:
#' \enumerate{
#'   \item An "inner integral": The part in the parentheses, which integrates over \eqn{p - 1}
#'   dimensions that are not part of the candidate splitting dimension \eqn{\ell = 1}. The result
#'   of this inner integral is a functioin of the one remaining variable, \eqn{x_1}. Call this
#'   resulting function \eqn{h(x_1)}.
#'   \item An "outer integral": A simple 1D integral of \eqn{h(x_1)} with respect to the splitting
#'   variable \eqn{x_1} over the domain of valid splitting thresholds \eqn{t \in [0,1]}.
#' }
#' The `precompute_integrals_along_dim` function is designed to numerically calculate this
#' inner integral function, \eqn{h(x_1)}, over a grid of points along the \eqn{x_1}-axis
#' (or, in general, the inner integral function \eqn{h(x_\ell)} over a grid of points along
#' the \eqn{x_\ell}-axis).
#'
#' @param ell The candidate splitting dimension along which we wish to pre-compute values of the
#'    inner integral function. See the **Details** section below.
#' @param p The total dimensionality of the variable of integration \eqn{x}.
#' @param theta_FUN_list The list of the \eqn{\theta^*(x) = (\theta_1^*(x),\ldots,\theta_K^*(x))}
#'    component functions \eqn{\theta_k^*:[0,1]^p\to\mathbb R}.
#' @param grid.size The number of points for the lookup table grid.
#' @param ... Additional arguments for the integrator passed to [cubature::hcubature()].
#'
#' @return An object of class `precomputed_integral` for fast querying.
#'
#' @importFrom cubature hcubature
#'
#' @export
precompute_integrals_along_dim <- function(ell, p, theta_FUN_list, grid.size = 200, ...) {
  K <- length(theta_FUN_list)
  grid.points <- seq(0, 1, length.out = grid.size)

  # The inner integrand
  vectorized_theta_FUN_p_minus_1 <- function(x.mat.p.minus.1, g) {
    n.pts <- if (p == 1) 1 else ncol(x.mat.p.minus.1)
    x.mat.p <- matrix(0.0, nrow = p, ncol = n.pts)
    x.ell.val <- g # The value along `ell`-th dimension is held constant

    # Reconstruct the full p-dimensional points before evaluation
    if (p > 1) {
      indices <- 1:(p-1)
      x.mat.p[-ell, ] <- x.mat.p.minus.1[indices, , drop = FALSE]
      x.mat.p[ell, ]  <- x.ell.val
    } else {
      x.mat.p[1, ] <- x.ell.val
    }

    x.t <- t(x.mat.p)
    theta.out <- lapply(theta_FUN_list, function(f) f(x.t))
    do.call(rbind, theta.out)
  }

  # Iterate over each point `g` on the grid and calculate the inner integral
  grid.int <- vapply(grid.points, function(g) {
    if (p == 1) {
      # Base case: Inner integral over 0 dimensions is the inner integral function itself
      vectorized_theta_FUN_p_minus_1(matrix(numeric(0), nrow = 0, ncol = 0), g = g)
    } else {
      # General case: Call hcubature for the (p-1)-dimensional inner integral
      res <- cubature::hcubature(
        f = function(x) vectorized_theta_FUN_p_minus_1(x, g = g),
        lowerLimit = rep(0, p - 1),
        upperLimit = rep(1, p - 1),
        fDim = K,
        vectorInterface = TRUE,
        ...
      )
      res$integral
    }
  }, FUN.VALUE = numeric(K))

  integral.value <- t(grid.int)
  colnames(integral.value) <- names(theta_FUN_list)

  structure(
    list(
      grid.points = grid.points,
      integral.value = integral.value,
      K = K,
      p = p,
      ell = ell
    ),
    class = "precomputed_integral"
  )
}

#' Queries the pre-computed lookup table to get parameter function integrals
#' across a candidate split
#'
#' Uses fast 1D interpolation and integration on the pre-computed "inner" integral
#' object return by [precompute_integrals_along_dim()] to find the intergal of
#' \eqn{\theta^*(x)} over two child regions defined by a binary partition of
#' the covariate space \eqn{[0,1]^p} into two hyperrectangles. The partition is
#' defined by splitting parameters \eqn{(\ell, t)}, where \eqn{\ell \in \{1,\ldots, p\}}
#' denotes the splitting axis/dimension along which the precomputed inner function
#' values were solved, and splitting threshold \eqn{t\in (0,1)}.
#'
#' @param precomputed_obj The object returned by [precompute_integrals_along_dim()].
#' @param t The splitting threshold \eqn{t \in (0, 1)} along the candidate splitting
#'    axis/dimension corresponding to pre-computed integral object `obj`.
#'
#' @return A list containing the integral values for the left child
#'    (\eqn{x_{\ell} \leq t}) and the right child (\eqn{x_{\ell} > t}) that
#'    result from the split along dimension \eqn{\ell \in \{1,\ldots, p\}} at
#'    threshold \eqn{t \in (0,1)}.
#'
#' @importFrom stats approxfun integrate
#'
#' @keywords internal
query_precomputed_integrals <- function(precomputed_obj, t) {
  K <- precomputed_obj$K
  integrals <- list(C1 = numeric(K), C2 = numeric(K))

  # Calculate the integrals for each of the K theta components
  for (k in 1:K) {
    # Create a fast 1D linear interpolation function from the discrete points
    interp_FUN <- stats::approxfun(
      x = precomputed_obj$grid.points,
      y = precomputed_obj$integral.value[, k],
      rule = 2 # `rule = 2` allows extrapolation to endpoints 0 and 1
    )
    # Use fast 1D integration on the interpolated function
    integrals$C1[k] <- stats::integrate(interp_FUN, lower = 0, upper = t)$value
    integrals$C2[k] <- stats::integrate(interp_FUN, lower = t, upper = 1)$value
  }
  return(integrals)
}

#' Compute population-level conditional expectations for the local parameters
#' across a candidate split
#'
#' Calculates the conditional expectation of the target effect functions
#' \eqn{\theta^*(x)} over the two child regions created by a binary split.
#' The calculation assumes that the covariates \eqn{X} are **uniformly distributed**
#' over the unit hypercube \eqn{[0,1]^p}. Under this assumption, the conditional
#' expectation simplifies to the integral of \eqn{\theta^*(x)} over the region,
#' divided by the region's volume.
#'
#' @details
#' This function computes \eqn{\theta_{C_j}^* := \mathbb{E}[\theta^*(X) \mid X \in S_j]}
#' for the two child regions \eqn{S_1, S_2} given split parameters \eqn{\ell, t}. The
#' split is defined by a dimension \eqn{\ell} and a threshold \eqn{t}, where
#' \deqn{S_1 \equiv S_1(\ell, t) = \{x \in [0,1]^p : x_\ell \leq t\}}
#' and
#' \deqn{S_2 \equiv S_2(\ell, t) = \{x \in [0,1]^p : x_\ell > t\}.}
#'
#' @inheritParams query_precomputed_integrals
#' @inheritParams precompute_integrals_along_dim
#' @param ... Additional arguments passed to the numerical integrator.
#'    See [precompute_integrals_along_dim()].
#'
#' @return A list containing two named numeric vectors, `C1` and `C2`,
#'   representing the conditional expectations \eqn{\theta_{C_1}^*} and
#'   \eqn{\theta_{C_2}^*} solved over the corresponding child region.
#'
#' @export
compute_theta_star <- function(precomputed_obj, t, theta_FUN_list, ...) {
  #--- Query the object to get the integral values over the two child regions
  integrals <- query_precomputed_integrals(precomputed_obj, t = t)

  #--- Calculate the conditional expectations by dividing by region volume.
  # For a split on dimension `ell` at `t` on the unit hypercube [0,1]^p:
  #   * Volume of C1(x_ell <= t) is t * 1^(p-1) = t
  #   * Volume of C2(x_ell > t) is (1-t) * 1^(p-1) = 1 - t
  vol.C1 <- t
  vol.C2 <- 1 - t

  K <- length(theta_FUN_list)
  theta.star.C1 <- numeric(K)
  theta.star.C2 <- numeric(K)

  # Avoid division by zero if t is at the boundary
  if (vol.C1 > sqrt(.Machine$double.eps)) theta.star.C1 <- integrals$C1 / vol.C1
  if (vol.C2 > sqrt(.Machine$double.eps)) theta.star.C2 <- integrals$C2 / vol.C2

  names(theta.star.C1) <- names(theta_FUN_list)
  names(theta.star.C2) <- names(theta_FUN_list)

  list(
    C1 = theta.star.C1,
    C2 = theta.star.C2
  )
}


