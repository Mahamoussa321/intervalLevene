#' Simulate interval data for k groups using a Gaussian copula
#'
#' Generates interval-valued data where each observation is defined by a center
#' and a radius. Dependence between center and radius is induced via a Gaussian
#' copula; margins can be Gamma or Poisson. For Poisson radii, zeros are bumped
#' to a tiny positive value to keep radii strictly > 0.
#'
#' @param group_sizes Integer vector. Number of observations per group.
#' @param dist Character string. Margin family: "gamma" (default) or "poisson".
#' @param center_par List of parameters for the center margin.
#'   - If \code{dist="gamma"}: \code{list(mean, scale)}. (Legacy: \code{disp} as alias of \code{scale}.)
#'   - If \code{dist="poisson"}: \code{list(mean)} used as \code{lambda}.
#' @param radius_par Same structure as \code{center_par} for the radius margin.
#' @param rho Numeric in (-1,1). Latent Gaussian correlation between center and radius.
#' @return A list of data frames (one per group) with columns: \code{center}, \code{radius}, \code{group}.
#' @examples
#' ## Example 1: Gamma margins with positive dependence
#' set.seed(2)
#' sim <- sim_data(
#'   group_sizes = c(10, 11, 9),
#'   dist        = "gamma",
#'   center_par  = list(mean = 5,  scale = 2),
#'   radius_par  = list(mean = 10, scale = 1),
#'   rho         = 0.5
#' )
#' head(do.call(rbind, sim))
#'
#' ## Example 2: Poisson margins with mild negative dependence
#' set.seed(123)
#' simP <- sim_data(
#'   group_sizes = c(6, 5, 7),
#'   dist        = "poisson",
#'   center_par  = list(mean = 5),   # lambda for centers
#'   radius_par  = list(mean = 8),   # lambda for radii
#'   rho         = -0.3
#' )
#' head(do.call(rbind, simP))
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats pnorm qgamma qpois
#' @export
sim_data <- function(group_sizes, dist = c("gamma","poisson"),
                     center_par = list(mean = 5,  scale = 2),
                     radius_par = list(mean = 10, scale = 1),
                     rho = 0.5) {
  dist <- match.arg(dist)
  stopifnot(is.numeric(group_sizes), length(group_sizes) >= 1L, all(group_sizes >= 1))
  stopifnot(is.numeric(rho), length(rho) == 1L, abs(rho) < 1)

  # helpers (with back-compat for 'disp')
  get_mean  <- function(par, nm) {
    if (!is.null(par$mean)) return(par$mean)
    stop(sprintf("Provide %s$mean", nm), call. = FALSE)
  }
  get_scale <- function(par, nm) {
    if (!is.null(par$scale)) return(par$scale)
    if (!is.null(par$disp)) {
      warning(sprintf("%s$disp is deprecated; use %s$scale instead.", nm, nm), call. = FALSE)
      return(par$disp)
    }
    stop(sprintf("Provide %s$scale (or legacy %s$disp) for Gamma margins.", nm, nm), call. = FALSE)
  }

  # parameter checks
  if (dist == "gamma") {
    if (get_mean(center_par, "center_par") <= 0 || get_scale(center_par, "center_par") <= 0)
      stop("Gamma center_par must have mean>0 and scale>0.", call. = FALSE)
    if (get_mean(radius_par, "radius_par") <= 0 || get_scale(radius_par, "radius_par") <= 0)
      stop("Gamma radius_par must have mean>0 and scale>0.", call. = FALSE)
  } else {
    if (get_mean(center_par, "center_par") < 0 || get_mean(radius_par, "radius_par") < 0)
      stop("Poisson means must be nonnegative.", call. = FALSE)
  }

  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)

  lapply(seq_along(group_sizes), function(i) {
    Z <- mvtnorm::rmvnorm(group_sizes[i], sigma = Sigma)  # (Z1,Z2) ~ N(0, Sigma)
    U <- stats::pnorm(Z)                                  # map to U(0,1)

    if (dist == "gamma") {
      thC <- get_scale(center_par, "center_par")
      thR <- get_scale(radius_par, "radius_par")
      kC  <- get_mean(center_par, "center_par") / thC
      kR  <- get_mean(radius_par, "radius_par") / thR
      xC  <- stats::qgamma(U[,1], shape = kC, scale = thC)
      xR  <- stats::qgamma(U[,2], shape = kR, scale = thR)
    } else {
      xC  <- stats::qpois(U[,1], lambda = get_mean(center_par, "center_par"))
      xR  <- stats::qpois(U[,2], lambda = get_mean(radius_par, "radius_par"))
      xR  <- pmax(xR, .Machine$double.eps)  # keep radius strictly > 0
    }

    data.frame(center = xC, radius = xR, group = i)
  })
}
