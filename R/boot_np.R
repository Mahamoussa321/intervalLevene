#' Nonparametric (case-resampling) bootstrap for the interval-Levene F
#'
#' Uses a **within-group nonparametric bootstrap**: resample rows *within each
#' group* (preserving group sizes) and recompute the interval-Levene F to get
#' an empirical reference for critical values and p-values.
#'
#' @param df List of per-group data.frames (first two cols center/radius),
#'   or one data.frame with `center`, `radius`, `group`.
#' @param omega Nonnegative weight on the radius component in Y. Default 1.
#' @param R Number of bootstrap replicates. Default 2000.
#' @param seed Optional integer seed.
#' @param alpha Significance level for the critical value (default 0.05).
#' @param return_boot If TRUE, return the full `boot::boot` object.
#'
#' @return list(F_obs, crit, p_value, R, boot?)
#' @export
#' @importFrom stats quantile
bootstrap_interval_levene_np <- function(df,
                                         omega = 1,
                                         R = 2000,
                                         seed = NULL,
                                         alpha = 0.05,
                                         return_boot = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  df <- .as_interval_df(df)

  stat_F <- function(d, i) {
    di <- d[i, , drop = FALSE]
    di$group <- factor(di$group, levels = levels(d$group))
    levene_interval_F(di, omega = omega)
  }

  F_obs <- stat_F(df, seq_len(nrow(df)))

  b <- boot::boot(
    data      = df,
    statistic = function(dat, ind) stat_F(dat, ind),
    R         = R,
    strata    = df$group
  )

  crit <- stats::quantile(b$t[, 1], 1 - alpha, names = FALSE)
  pval <- mean(b$t[, 1] >= F_obs)

  out <- list(F_obs = F_obs, crit = crit, p_value = pval, R = R)
  if (return_boot) out$boot <- b
  out
}
