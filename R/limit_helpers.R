#' Bind a list of per-group data.frames into one long data.frame
#'
#' Each element of `sim_list` is a group's data with columns `center` and `radius`.
#' This function adds a `group` factor and row-binds them.
#'
#' @param sim_list list of data.frames (one per group) with columns `center`, `radius`
#' @return data.frame with columns `center`, `radius`, `group` (factor)
#' @export
#' @importFrom data.table rbindlist
il_bind_groups <- function(sim_list) {
  stopifnot(is.list(sim_list), length(sim_list) >= 1L)
  df <- data.table::rbindlist(
    lapply(seq_along(sim_list), function(i) {
      x <- sim_list[[i]]
      if (!all(c("center","radius") %in% names(x))) {
        stop("Each element must have columns 'center' and 'radius'.", call. = FALSE)
      }
      transform(x, group = i)
    }),
    use.names = TRUE, fill = TRUE
  )
  df$group <- factor(df$group)
  df
}

#' Pooled (within-group) moments for centers/radii
#'
#' Computes pooled within-group variance of `center` and `radius`, and their
#' pooled within-group correlation using residuals from each group's mean.
#'
#' @param df data.frame with columns `center`, `radius`, `group` (factor)
#' @return list(s2C = var(center residuals), s2R = var(radius residuals), rho = cor residuals)
#' @export
#' @importFrom stats cor
il_pool_moments <- function(df) {
  stopifnot(is.data.frame(df), all(c("center","radius","group") %in% names(df)))
  df$group <- as.factor(df$group)
  muC_g <- tapply(df$center, df$group, mean)
  muR_g <- tapply(df$radius, df$group, mean)
  resC  <- df$center - muC_g[as.character(df$group)]
  resR  <- df$radius - muR_g[as.character(df$group)]
  s2C <- mean(resC^2)
  s2R <- mean(resR^2)
  rho <- stats::cor(resC, resR)
  list(s2C = as.numeric(s2C), s2R = as.numeric(s2R), rho = as.numeric(rho))
}

#' Limiting draws for the interval-Levene F under H0
#'
#' Simulates the asymptotic null reference:
#' \deqn{ F_\infty = \alpha \,\frac{\chi^2_{df,C}}{df} + \beta \,\frac{\chi^2_{df,R}}{df}, }
#' where \eqn{\alpha = \sigma_C^2/(\sigma_C^2 + \omega \sigma_R^2)},
#' \eqn{\beta = \omega \sigma_R^2/(\sigma_C^2 + \omega \sigma_R^2)},
#' and the two chi-squares arise from sums of squares of correlated normals with correlation \eqn{\rho}.
#'
#' @param k integer, number of groups
#' @param M integer, number of draws
#' @param omega nonnegative numeric weight
#' @param s2C pooled within-group variance of centers
#' @param s2R pooled within-group variance of radii
#' @param rho pooled within-group correlation between center and radius residuals
#' @return numeric vector of length `M` with limiting draws
#' @export
#' @importFrom mvtnorm rmvnorm
il_limit_draws <- function(k, M, omega, s2C, s2R, rho) {
  stopifnot(is.numeric(k), k >= 2, is.numeric(M), M >= 1,
            is.numeric(omega), omega >= 0)
  df <- k - 1L
  denom <- s2C + omega * s2R
  if (!(is.finite(denom) && denom > 0)) {
    stop("Nonpositive or non-finite denominator (s2C + omega*s2R).", call. = FALSE)
  }
  alpha <- s2C / denom
  beta  <- (omega * s2R) / denom

  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  out <- numeric(M)
  for (m in seq_len(M)) {
    Z  <- mvtnorm::rmvnorm(df, sigma = Sigma)  # df x 2
    SC <- sum(Z[, 1]^2)
    SR <- sum(Z[, 2]^2)
    out[m] <- alpha * (SC / df) + beta * (SR / df)
  }
  out
}
