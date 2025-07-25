#' Create a criterion definition object
#'
#' This constructor function creates a standardized object that defines a splitting
#' criterion. This allows the main simulation loop to be agnostic to the specific
#' calculations of any given criterion.
#'
#' @param name A character string, the name of the criterion.
#' @param precompute A function that takes the full simulation data list (`sim.data`)
#'    and performs any necessary calculations on the parent node before the split
#'    loop begins. It should return a list or object containing pre-computed values.
#'    If no pre-computation is needed, this can be `NULL`.
#' @param compute A function that calculates the criterion's value for a single split.
#'    It must accept two arguments:
#'    \describe{
#'      \item{`split.data`}{A list containing data for the current split, such as
#'        child node estimates (`theta.hat.C1`, `theta.hat.C2`), population
#'        parameters (`theta.star.C1`, `theta.star.C2`), the weight (`wt`), and the
#'        Jacobian/covariance matrix (`V`).}
#'     \item{`precomputed.data`}{The object returned by the `precompute` function.}
#'   }
#'
#' @return An object of class `criteval_criterion`, which is a list containing the
#'   `name`, `precompute`, and `compute` functions.
#'
#' @keywords internal
new_criterion <- function(name, precompute = NULL, compute) {
  stopifnot(is.character(name), is.function(compute))
  if (!is.null(precompute)) stopifnot(is.function(precompute))

  structure(
    list(
      name = name,
      precompute = precompute,
      compute = compute
    ),
    class = "criteval_criterion"
  )
}

#' Get a list of predefined criterion objects
#'
#' This function returns a named list of predefined `criteval_criterion` objects.
#' This serves as a central registry for all available criteria, making it easy
#' to extend the framework by simply adding new definitions here.
#'
#' @param names An optional character vector of criterion names to retrieve.
#'   If `NULL` (the default), a named list of all available criteria is returned.
#'
#' @return A named list of `criteval_criterion` objects. If `names` is a single
#'   string, a single `criteval_criterion` object is returned.
#'
#' @examples
#' \dontrun{
#' get_criteria()          # Returns the full list
#' get_criteria("Delta.V") # Returns just the Delta.V object
#' get_criteria(c("Delta", "Delta.pop")) # Returns a list with just those two
#' get_criteria("foo")     # Throws a helpful error
#' }
#'
#' @export
get_criteria <- function(names = NULL) {
  .registry <- list(
    new_criterion(
      name = "Delta.pop",
      compute = function(split.data, precomputed.data) {
        diff.star <- split.data$theta.star.C1 - split.data$theta.star.C2
        split.data$wt * sum(diff.star^2)
      }
    ),
    new_criterion(
      name = "Delta.V.pop",
      compute = function(split.data, precomputed.data) {
        diff.star <- split.data$theta.star.C1 - split.data$theta.star.C2
        split.data$wt * sum((split.data$V %*% diff.star)^2)
      }
    ),
    new_criterion(
      name = "Delta",
      compute = function(split.data, precomputed.data) {
        diff.hat <- split.data$theta.hat.C1 - split.data$theta.hat.C2
        split.data$wt * sum(diff.hat^2)
      }
    ),
    new_criterion(
      name = "Delta.V",
      compute = function(split.data, precomputed.data) {
        diff.hat <- split.data$theta.hat.C1 - split.data$theta.hat.C2
        split.data$wt * sum((split.data$V %*% diff.hat)^2)
      }
    ),
    new_criterion(
      name = "Delta.grad",
      precompute = function(sim.data) {
        W <- scale(sim.data$W, scale = FALSE)
        Y <- scale(sim.data$Y, scale = FALSE)
        n <- nrow(W)

        # theta.hat.P <- solve(crossprod(W)) %*% crossprod(W, Y)
        theta.hat.P <- solve_OLS(w.mat = W, y.vec = Y, intercept = F)
        psi.theta <- sweep(W, 1, Y - W %*% theta.hat.P, "*")
        A.P.inv <- solve(crossprod(W)/n)

        list(rho.grad = psi.theta %*% A.P.inv)
      },
      compute = function(split.data, precomputed.data) {
        rho.grad <- precomputed.data$rho.grad

        diff.rho.grad <-
          colMeans(rho.grad[split.data$C1, , drop = FALSE]) -
          colMeans(rho.grad[split.data$C2, , drop = FALSE])

        split.data$wt * sum(diff.rho.grad^2)
      }
    ),
    new_criterion(
      name = "Delta.fpt",
      precompute = function(sim.data) {
        W <- scale(sim.data$W, scale = FALSE)
        Y <- scale(sim.data$Y, scale = FALSE)

        # theta.hat.P <- solve(crossprod(W)) %*% crossprod(W, Y)
        theta.hat.P <- solve_OLS(w.mat = W, y.vec = Y, intercept = F)
        psi.theta <- sweep(W, 1, Y - W %*% theta.hat.P, "*")

        list(rho.fpt = psi.theta)
      },
      compute = function(split.data, precomputed.data) {
        rho.fpt <- precomputed.data$rho.fpt

        diff.rho.fpt <-
          colMeans(rho.fpt[split.data$C1, , drop = FALSE]) -
          colMeans(rho.fpt[split.data$C2, , drop = FALSE])

        split.data$wt * sum(diff.rho.fpt^2)
      }
    )
  )

  .registry <- stats::setNames(.registry, nm = sapply(.registry, function(x) x$name))
  if (is.null(names)) {
    return(.registry)
  }

  unmatched <- setdiff(names, names(.registry))
  if (length(unmatched) > 0L) {
    stop("The following criteria were not found: ", paste(unmatched, collapse = ", "))
  }

  result <- .registry[names]
  if (length(result) == 1L) {
    return(result[[1L]])
  }

  return(result)
}
