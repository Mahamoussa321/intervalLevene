#' Interval-Levene F statistic
#'
#' Computes an ANOVA-style Levene statistic on interval-valued data by forming
#' \eqn{Y_{ij} = (C_{ij}-\bar C_i)^2 + \omega (R_{ij}-\bar R_i)^2} and returning
#' \eqn{F = \{SSB/(k-1)\}/\{SSE/(N-k)\}}.
#'
#' @param df data.frame with columns \code{center}, \code{radius}, and (optionally) \code{group}.
#' @param group Factor (or coercible) of length \code{nrow(df)} giving group memberships.
#'   If \code{NULL} and \code{df$group} exists, it will be used.
#' @param omega Nonnegative numeric weight on the radius component (default 1).
#' @return Numeric F statistic (MSB/MSE). Returns \code{NA_real_} with a warning if SSE=0.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   center = rnorm(60),
#'   radius = rexp(60),
#'   group  = factor(rep(1:3, each = 20))
#' )
#' levene_interval_F(df, df$group)
#' @export
levene_interval_F <- function(df, group = NULL, omega = 1) {
  # basic checks
  stopifnot(is.data.frame(df))
  if (!all(c("center", "radius") %in% names(df)))
    stop("df must contain columns 'center' and 'radius'.", call. = FALSE)

  if (is.null(group)) {
    if ("group" %in% names(df)) {
      group <- df$group
    } else {
      stop("Provide 'group' or include a 'group' column in df.", call. = FALSE)
    }
  }
  df$group <- as.factor(group)
  k <- nlevels(df$group); N <- nrow(df)
  if (k < 2L) stop("Need at least 2 groups.", call. = FALSE)
  if (!is.numeric(omega) || length(omega) != 1L || omega < 0)
    stop("'omega' must be a single nonnegative number.", call. = FALSE)

  # drop NAs (or choose to error instead)
  rows_ok <- stats::complete.cases(df$center, df$radius, df$group)
  if (!all(rows_ok)) {
    df <- df[rows_ok, , drop = FALSE]
    N  <- nrow(df)
    warning("Missing values removed before computing the statistic.", call. = FALSE)
  }

  # per-group means of center and radius
  muC_g <- dplyr::summarise(dplyr::group_by(df, group),
                            muC = mean(center), .groups = "drop")
  muR_g <- dplyr::summarise(dplyr::group_by(df, group),
                            muR = mean(radius), .groups = "drop")
  df_g  <- dplyr::left_join(df, muC_g, by = "group") |>
    dplyr::left_join(muR_g, by = "group")

  # transformed response
  Y   <- (df_g$center - df_g$muC)^2 + omega * (df_g$radius - df_g$muR)^2
  dfY <- data.frame(group = df_g$group, Y = Y)

  # group means of Y and sizes
  Ybar_g <- dplyr::summarise(dplyr::group_by(dfY, group),
                             Ybar = mean(Y), n = dplyr::n(), .groups = "drop")
  if (any(Ybar_g$n < 2L))
    warning("Some groups have size 1; resulting F may be unstable.", call. = FALSE)

  Ybar <- mean(dfY$Y)

  # sums of squares
  SSB <- sum(Ybar_g$n * (Ybar_g$Ybar - Ybar)^2)   # between
  SSE <- dfY |>
    dplyr::left_join(Ybar_g[, c("group","Ybar")], by = "group") |>
    dplyr::summarise(SSE = sum((Y - Ybar)^2)) |>
    dplyr::pull(SSE)

  df1 <- k - 1
  df2 <- N - k
  if (df2 <= 0) stop("Not enough total observations: need N > k.", call. = FALSE)
  if (SSE <= 0) {
    warning("SSE is zero; returning NA.", call. = FALSE)
    return(NA_real_)
  }

  # F = MSB/MSE
  (SSB / df1) / (SSE / df2)
}
