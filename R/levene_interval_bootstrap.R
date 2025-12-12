#' Bootstrap calibration for interval-valued Levene-type tests
#'
#' Implements Algorithm 1 (Ramos–Lubiano, 2012) for interval data.
#' Works on either a data.frame with columns `L`, `R`, and a grouping
#' column, or on a list of per-group data.frames with columns `L` and `R`.
#'
#' @param data A data.frame (with columns L,R and a grouping column) or
#'   a list of data.frames (one per group) each containing columns L and R.
#' @param group For data.frames only: the name (string) of the grouping column.
#' @param theta Nonnegative weight \eqn{\omega} for the radius term. Default 1.
#' @param B Number of bootstrap resamples. Default 2000.
#' @param version Either "T1" or "T2". See paper for definitions.
#' @param alpha Test size for the decision rule \eqn{q_{1-\alpha}}. Default 0.05.
#' @param seed Optional integer seed (reproducible resampling). Default NULL.
#' @param return_draws If TRUE, include the full vector of bootstrap draws.
#' @return A list with elements:
#'   \item{method}{Character, "bootstrap T1" or "bootstrap T2".}
#'   \item{stat}{Observed statistic \eqn{T_{\text{obs}}}.}
#'   \item{q_alpha}{Bootstrap critical value \eqn{q_{1-\alpha}}.}
#'   \item{p.value}{Bootstrap p-value \eqn{\Pr(T^\ast \ge T_{\text{obs}})}.}
#'   \item{reject}{Logical; TRUE if \eqn{T_{\text{obs}} > q_{1-\alpha}}.}
#'   \item{B, theta, version, alpha}{Echoed inputs.}
#'   \item{draws}{(optional) length-\code{B} vector of bootstrap draws.}
#' @export
#' @examples
#' \dontrun{
#' # df: columns L,R,group
#' out <- levene_interval_bootstrap(df, group = "group", theta = 1, B = 2000, version = "T1")
#' out$reject; out$p.value
#' }
levene_interval_bootstrap <- function(
    data, group = NULL, theta = 1, B = 2000,
    version = c("T1", "T2"), alpha = 0.05,
    seed = NULL, return_draws = FALSE
){
  version <- match.arg(version)

  # -- Normalize to list-of-groups with L,R columns --------------------------
  X <- if (is.list(data) && !is.data.frame(data)) {
    # list-of-dfs path
    lapply(data, function(d) {
      stopifnot(is.data.frame(d), all(c("L","R") %in% names(d)))
      d <- transform(d, L = as.numeric(L), R = as.numeric(R))
      swap <- which(d$L > d$R); if (length(swap)) { tmp <- d$L[swap]; d$L[swap] <- d$R[swap]; d$R[swap] <- tmp }
      d
    })
  } else {
    # data.frame path
    stopifnot(is.data.frame(data))
    if (is.null(group)) stop("For data.frames, supply `group` = name of grouping column.")
    if (!all(c("L","R", group) %in% names(data))) {
      stop("`data` must have columns L, R, and the specified `group`.")
    }
    split(transform(data, L = as.numeric(L), R = as.numeric(R)), data[[group]])
  }

  k   <- length(X)
  n_i <- vapply(X, nrow, integer(1))
  stopifnot(k >= 2L, all(n_i >= 2L))
  if (!is.null(seed)) set.seed(seed)
  N <- sum(n_i)

  # -- Precompute from original data -----------------------------------------
  Ybar_i <- numeric(k)   # \bar Y_{i·}
  den_i  <- numeric(k)   # (1/n_i) Σ (Y_ij - \bar Y_{i·})^2
  mu_m   <- mu_s <- numeric(k)

  for (i in seq_len(k)) {
    d <- X[[i]]
    m <- (d$L + d$R)/2
    s <- (d$R - d$L)/2
    mbar <- mean(m); sbar <- mean(s)
    mu_m[i] <- mbar;  mu_s[i] <- sbar
    Y <- (m - mbar)^2 + theta*(s - sbar)^2
    Ybar_i[i] <- mean(Y)
    den_i[i]  <- mean((Y - mean(Y))^2)
  }
  Ydotdot  <- sum(n_i * Ybar_i) / N
  T_obs    <- sum(n_i * (Ybar_i - Ydotdot)^2) / sum(den_i)
  T1_den   <- sum(den_i)

  # -- Bootstrap loop ---------------------------------------------------------
  draws <- numeric(B)
  for (b in seq_len(B)) {
    Ybar_star_i <- numeric(k)
    den_star_T2 <- 0

    for (i in seq_len(k)) {
      d <- X[[i]]
      idx <- sample.int(n_i[i], n_i[i], replace = TRUE)
      Lst <- d$L[idx]; Rst <- d$R[idx]
      swap <- which(Lst > Rst); if (length(swap)) { tmp <- Lst[swap]; Lst[swap] <- Rst[swap]; Rst[swap] <- tmp }
      m_st <- (Lst + Rst)/2
      s_st <- (Rst - Lst)/2

      Yst <- (m_st - mu_m[i])^2 + theta*(s_st - mu_s[i])^2
      Ybar_star_i[i] <- mean(Yst)

      if (identical(version, "T2")) {
        den_star_T2 <- den_star_T2 +
          (1 / (n_i[i] * pmax(Ybar_i[i]^2, .Machine$double.eps))) *
          sum((Yst - mean(Yst))^2)
      }
    }

    delta_i <- Ybar_star_i - Ybar_i
    if (identical(version, "T1")) {
      w_i   <- Ydotdot / pmax(Ybar_i, .Machine$double.eps)
      grand <- (1 / N) * sum(n_i * w_i * delta_i)
      num   <- sum( n_i * (w_i * delta_i - grand)^2 )
      denom <- T1_den
    } else {
      w_i   <- 1 / pmax(Ybar_i, .Machine$double.eps)
      grand <- (1 / N) * sum(n_i * w_i * delta_i)
      num   <- sum( n_i * (w_i * delta_i - grand)^2 )
      denom <- den_star_T2
    }
    draws[b] <- num / denom
  }

  q_alpha <- as.numeric(stats::quantile(draws, 1 - alpha, names = FALSE, type = 7))
  pval    <- mean(draws >= T_obs)
  out <- list(
    method  = paste("bootstrap", version),
    stat    = T_obs,
    q_alpha = q_alpha,
    p.value = pval,
    reject  = (T_obs > q_alpha),
    B = B, theta = theta, version = version, alpha = alpha
  )
  if (isTRUE(return_draws)) out$draws <- draws
  out
}
