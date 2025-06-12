#' Run a customizable splitting criteria simulation under a varying-coefficient
#' model design
#'
#' Orchestrates a simulation to evaluate and compare splitting criteria in their
#' ability to detect heterogeneity in an unobserved parameter \eqn{\theta^*(x)}.
#' The task of the criteria is to select a single binary axis-aligned split of
#' the data such as to separate the covariate space into dissimilar regions
#' with respect to \eqn{\theta^*(x)}.
#'
#' @param n Integer, the total number of observations.
#' @param kappa.ratio Numeric, the target condition number ratio for the `W` covariance matrix.
#' @param rho.W Numeric, the AR(1) correlation for `W` regressors.
#' @param R.W A pre-specified correlation matrix for `W`.
#' @param q.W Integer, the number of "strong" (high-variance) `W` regressors.
#' @param mu.W Numeric vector, the mean of `W` regressors.
#' @param p.X Integer, the dimension of the auxiliary covariates `X`.
#' @param sigma.eps Numeric, the standard deviation of the error term.
#' @param theta_FUN_list A list of functions defining the `theta*(x)` components.
#' @param nu_FUN A function defining the nuisance component `nu(x)`.
#' @param criteria A list of `criteval_criterion` objects to evaluate. Defaults to all
#'    criteria from `criteval::get_criteria()`.
#' @param min.node.size Integer, the minimum number of observations in any child node.
#' @param .n.skip Integer, interval for skipping split points to speed up computation.
#'    Computes criteria every `.n.skip` values. Defaults to 1.
#' @param .grid.size Integer, the resolution for the numerical integration lookup table.
#' @param .tol.integrate Tolerance parameter passed to [cubature::hcubature()] when
#'    performing numerical integration.
#' @param seed Integer, an optional random seed for reproducibility.
#'
#' @note The numerical integration for \eqn{theta^*(x)} assumes that the auxiliary
#'   covariates \eqn{X} are uniformly distributed over the covariate space \eqn{[0,1]^p}.
#'   Therefore, the `rho.X` parameter for the `data_generator()` is intentionally not
#'   exposed because the general Gaussian copula would violate these assumptions and
#'   lead to spurrious expectations \eqn{\mathbb E[\theta^*(X) \mid X \in C_j]}.
#'
#' @export
#'
#' @importFrom dplyr %>% .data bind_rows mutate group_by ungroup filter recode
run_criteval_sim <- function(
    n,
    kappa.ratio, rho.W = NULL, R.W = NULL, q.W, mu.W = NULL,
    p.X, sigma.eps,
    theta_FUN_list, nu_FUN,
    criteria = criteval::get_criteria(),
    min.node.size = 10,
    .n.skip = 1, # Only compute criterion values every .n.skip values
    .grid.size = 200,
    .tol.integrate = 1e-4,
    seed = NULL) {
  #--- Validation
  stopifnot(is.list(criteria) && all(sapply(criteria, inherits, "criteval_criterion")))
  crit.names <- sapply(criteria, `[[`, "name")

  #--- Generation
  message(sprintf("%s :: %s...", format(Sys.time()), "Generating data"))
  dgp <- criteval::data_generator(
    kappa.ratio    = kappa.ratio,
    rho.W          = rho.W,
    R.W            = R.W,
    q.W            = q.W,
    mu.W           = mu.W,
    p.X            = p.X,
    sigma.eps      = sigma.eps,
    theta_FUN_list = theta_FUN_list,
    nu_FUN         = nu_FUN
  )
  sim.data <- dgp$generate(n, seed = seed)

  #--- "Parent-level" pre-computation phase. Calculations done once over the initial set
  message(sprintf("%s :: %s...", format(Sys.time()), "Parent-level computation phase for criteria"))
  precomputed.data <- lapply(criteria, function(crit) {
    if (!is.null(crit$precompute)) {
      crit$precompute(sim.data)
    } else {
      NULL
    }
  })

  #--- Primary loop over the sequence of candidate splits
  # TODO: Make sure this is always an increasing sequence, non-generate, non-empty, etc.
  split.indices <- seq(min.node.size, n - min.node.size, by = .n.skip)
  res.by.dim <- vector(mode = "list", length = p.X)

  # Outer loop over covariate dimensions
  for (split.dim in seq_len(p.X)) {
    message(sprintf("%s :: Scanning covariate feature %*i of %i...", format(Sys.time()), nchar(p.X), split.dim, p.X))

    precomputed_integrals <- precompute_integrals_along_dim(
      ell = split.dim,
      p = p.X,
      theta_FUN_list = theta_FUN_list,
      grid.size = .grid.size,
      tol = .tol.integrate
    )

    obs.order <- order(sim.data$X[, split.dim])

    # Inner loop over candidate thresholds along the dimension `split.dim`
    res.by.thresh <- sapply(split.indices, function(i) {
      split.thresh <- sim.data$X[obs.order[i], split.dim]
      C1 <- sim.data$X[, split.dim] <= split.thresh
      C2 <- !C1

      n.C1 <- sum(C1)
      n.C2 <- sum(C2)

      if (n.C1 < min.node.size || n.C2 < min.node.size) {
        return(structure(NaN, length(criteria)), names = crit.names)
      }

      # Population parameters/conditional expectation over the child
      theta.star <- compute_theta_star(
        precomputed_obj = precomputed_integrals,
        t = split.thresh,
        theta_FUN_list = theta_FUN_list
      )
      # Empirical solution/finite-sample estimates
      theta.hat <- compute_theta_hat(
        C1 = C1,
        C2 = C2,
        w.mat = sim.data$W,
        y.vec = sim.data$Y,
        labels = names(theta_FUN_list)
      )

      split.data <- list(
        "wt" = n.C1 * n.C2 / (n.C1 + n.C2)^2,
        "C1" = which(C1),
        "C2" = which(C2),
        "V"  = dgp$Sigma.W,
        "theta.star.C1" = theta.star$C1,
        "theta.star.C2" = theta.star$C2,
        "theta.hat.C1" = theta.hat$C1,
        "theta.hat.C2" = theta.hat$C2
      )

      # Compute criteria values for this split
      sapply(criteria, function(crit) {
        crit$compute(split.data, precomputed.data[[crit$name]])
      })
    })

    res.by.dim[[split.dim]] <- data.frame(
      index = obs.order[split.indices],
      dim = split.dim,
      t(res.by.thresh)
    )
  }

  #--- Post-processing
  message(sprintf("%s :: Post-processing results...", format(Sys.time())))

  df.res <- res.by.dim %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(unname(crit.names)),
      names_to = "criterion",
      values_to = "value"
    ) %>%
    dplyr::mutate(criterion = factor(.data$criterion, levels = crit.names)) %>%
    dplyr::mutate(threshold = sim.data$X[cbind(.data$index, .data$dim)]) %>%
    dplyr::group_by(.data$criterion, .data$dim) %>%
    dplyr::mutate(
      value_norm = (.data$value - min(.data$value, na.rm = T))/(max(.data$value, na.rm = T) - min(.data$value, na.rm = T))
    ) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(
      cols = c("value", "value_norm"),
      names_to = "normalized",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      normalized = dplyr::recode(.data$normalized, value = "Raw", value_norm = "Normalized")
    ) %>%
    dplyr::mutate(
      normalized = factor(.data$normalized, levels = c("Normalized", "Raw"))
    ) %>%
    dplyr::filter(!is.nan(.data$value))
  return(df.res)
}
