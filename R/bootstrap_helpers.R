#' Residual bootstrap under H0 for the interval-Levene F
#'
#' Computes a Levene-type F on interval data where each observation is an
#' interval \eqn{[C \pm R]} via the proxy
#' \deqn{Y_{ij}=(C_{ij}-\bar C_i)^2 + \omega(R_{ij}-\bar R_i)^2,}
#' then calibrates the test using a **residual bootstrap under the null**
#' (equal group variances) by rescaling group residuals to the pooled SD
#' and resampling within groups.
#'
#' @param df A list of per-group data.frames (first two columns are center/radius),
#'   or a single data.frame with columns `center`, `radius`, `group` (factor-able).
#' @param omega Nonnegative weight on the radius term in \eqn{Y}. Default 1.
#' @param R Number of bootstrap replicates. Default 2000.
#' @param seed Optional integer seed for reproducibility.
#' @param return_draws If `TRUE`, include the vector of bootstrap draws `F_star`.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{F_obs}: observed interval-Levene F statistic.
#'   \item \code{crit_95}: 0.95 bootstrap critical value for \code{F}.
#'   \item \code{p_value}: bootstrap p-value, \eqn{\mathrm{mean}(F^\ast \ge F_{\mathrm{obs}})}.
#'   \item \code{R}: number of bootstrap replicates used.
#'   \item \code{F_star} (present if \code{return_draws=TRUE}): vector of bootstrap draws.
#' }
#'
#' @details
#' The bootstrap “imposes the null” of equal variances by (i) removing group means,
#' (ii) rescaling each group’s residuals to the pooled SD, (iii) adding back the
#' common mean, and (iv) resampling within groups to original sizes. This yields a
#' bootstrap null for the F statistic that respects group sizes and mean structure.
#'
#' @examples
#' \donttest{
#' ## Example 1: Single data.frame input (3 groups, unequal variances)
#' set.seed(123)
#' n_per_group <- c(40, 45, 50)
#' g <- factor(rep(LETTERS[1:3], times = n_per_group))
#' muC <- c(A = 0,  B = 0.5, C = -0.5)[g]
#' muR <- c(A = 1,  B = 1,   C = 1)[g]
#' sdC <- c(A = 0.5, B = 0.8, C = 1.1)[g]
#' sdR <- c(A = 0.4, B = 0.4, C = 0.9)[g]
#' dat <- data.frame(
#'   center = rnorm(length(g), mean = muC, sd = sdC),
#'   radius = pmax(1e-6, rnorm(length(g), mean = muR, sd = sdR)),
#'   group  = g
#' )
#' fit1 <- bootstrap_interval_levene(dat, omega = 1, R = 399, seed = 1)
#' fit1$F_obs; fit1$crit_95; fit1$p_value
#' # Reject H0 at 5% if F_obs > crit_95 (equivalently p_value < 0.05)
#'
#' ## Example 2: List-of-data.frames input (same data split by group)
#' lst <- split(dat[c("center","radius")], dat$group)
#' fit2 <- bootstrap_interval_levene(lst, omega = 0.5, R = 399, seed = 1)
#' fit2$p_value
#'
#' ## Example 3: Return bootstrap draws to visualize the null distribution
#' fit3 <- bootstrap_interval_levene(dat, omega = 1, R = 499, seed = 42, return_draws = TRUE)
#' hist(fit3$F_star, breaks = 30, main = "Bootstrap null of F", xlab = "F*")
#' abline(v = fit3$F_obs, lwd = 2)
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{car::leveneTest()} for classical Levene/Brown–Forsythe tests.
#'   \item \code{stats::var.test()} for two-sample F test on variances (parametric).
#' }
#'
#' @references
#' Levene, H. (1960). Robust tests for equality of variances. In \emph{Contributions
#' to Probability and Statistics: Essays in Honor of Harold Hotelling}, 278–292.
#'
#' Brown, M.B., & Forsythe, A.B. (1974). Robust tests for equality of variances.
#' \emph{Journal of the American Statistical Association}, 69(346), 364–367.
#'
#' @export
#' @importFrom stats sd quantile

bootstrap_interval_levene <- function(df, omega = 1, R = 2000, seed = NULL,
                                      return_draws = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  df <- .as_interval_df(df)             # center, radius, group

  dY  <- .build_Y(df, omega)
  k   <- nlevels(dY$group)
  N   <- nrow(dY)
  if (k < 2L) stop("Need at least 2 groups.", call. = FALSE)

  F_obs <- .levene_F_from_Y(dY)

  # --- residual bootstrap under H0 (equal Var(Y) across groups) ---
  Ybar_g <- tapply(dY$Y, dY$group, mean)
  n_g    <- as.numeric(table(dY$group))
  Ybar   <- mean(dY$Y)
  resid  <- dY$Y - rep(Ybar_g, times = n_g)

  s_i  <- tapply(resid, dY$group, stats::sd)           # group sd
  s2_p <- sum((n_g - 1) * s_i^2) / (N - k)             # pooled var
  s_p  <- sqrt(s2_p)
  scale_by <- s_p / s_i
  scale_by[!is.finite(scale_by)] <- 1                  # handle sd=0

  Y0 <- Ybar + resid * rep(scale_by, times = n_g)      # null-imposed Y
  Y0_split <- split(Y0, dY$group)
  g_levels <- levels(dY$group)

  F_star <- numeric(R)
  for (b in seq_len(R)) {
    # resample within groups to original sizes
    Yb_split <- mapply(function(v, m) sample(v, size = m, replace = TRUE),
                       Y0_split, n_g, SIMPLIFY = FALSE)
    Yb_all  <- unlist(Yb_split, use.names = FALSE)
    grp_all <- factor(rep(g_levels, times = n_g), levels = g_levels)
    dYb <- data.frame(group = grp_all, Y = Yb_all)

    F_star[b] <- .levene_F_from_Y(dYb)
  }

  crit_95 <- stats::quantile(F_star, 0.95, names = FALSE)
  pval    <- mean(F_star >= F_obs)

  out <- list(F_obs = F_obs, crit_95 = crit_95, p_value = pval, R = R)
  if (return_draws) out$F_star <- F_star
  out
}

# -----------------------
# internal helpers
# -----------------------

# Coerce to data.frame(center, radius, group)
#' @keywords internal
.as_interval_df <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    df <- do.call(rbind, lapply(seq_along(x), function(i) {
      xi <- x[[i]]
      nm <- names(xi)
      if (is.null(nm) || !("center" %in% nm && "radius" %in% nm)) {
        names(xi)[1:2] <- c("center","radius")
      }
      xi$group <- i
      xi
    }))
    df$group <- factor(df$group)
    df[, c("center","radius","group")]
  } else if (is.data.frame(x)) {
    df <- x
    if (!"center" %in% names(df)) names(df)[1] <- "center"
    if (!"radius" %in% names(df)) names(df)[2] <- "radius"
    if (!"group"  %in% names(df)) stop("Provide a 'group' column.")
    df$group <- factor(df$group)
    df[, c("center","radius","group")]
  } else {
    stop("Input must be a list of per-group data.frames or a single data.frame.")
  }
}

# Build Y = (C - Cbar_i)^2 + omega (R - Rbar_i)^2
#' @keywords internal
.build_Y <- function(df, omega = 1) {
  muC_g <- tapply(df$center, df$group, mean)
  muR_g <- tapply(df$radius, df$group, mean)
  rC <- df$center - muC_g[df$group]
  rR <- df$radius - muR_g[df$group]
  Y  <- rC^2 + omega * rR^2
  data.frame(group = df$group, Y = Y)
}

# Compute F = MSB/MSE from data.frame(group, Y)
#' @keywords internal
.levene_F_from_Y <- function(dY) {
  k  <- nlevels(dY$group); N <- nrow(dY)
  Ybar_g <- tapply(dY$Y, dY$group, mean)
  n_g    <- as.numeric(table(dY$group))
  Ybar   <- mean(dY$Y)

  SSB <- sum(n_g * (Ybar_g - Ybar)^2)
  SSE <- sum((dY$Y - rep(Ybar_g, times = n_g))^2)

  if (!is.finite(SSE) || SSE <= 0) return(NA_real_)
  (SSB / (k - 1)) / (SSE / (N - k))
}
