#' Bootstrap T for interval data (Ramos–Lubiano 2012)
#'
#' Exact-null bootstrap for Levene-type T (T1/T2) on Y = (m - mbar)^2 + omega (s - sbar)^2.
#' @param X list of data.frames (one per group) each with numeric columns L and R.
#' @param theta numeric scalar, weight on the radius component (use 1 for standard).
#' @param B integer, number of bootstrap replicates.
#' @param version character, "T1" or "T2".
#' @param seed optional integer seed for reproducibility (NULL to use ambient RNG).
#' @return list with elements \code{stat}, \code{p.value}, and \code{draws}.
#' @export
interval_bootstrap_RL2012 <- function(
    X, theta = 1, B = 2000, version = c("T1","T2"), seed = NULL
){
  version <- match.arg(version)
  if (!is.null(seed)) set.seed(seed)

  # ---- Input checks & normalization ----
  stopifnot(is.list(X), length(X) >= 2L)
  k    <- length(X)
  n_i  <- sapply(X, nrow)
  stopifnot(all(n_i >= 2))
  for (i in seq_len(k)) {
    stopifnot(all(c("L","R") %in% names(X[[i]])))
    bad <- which(X[[i]]$L > X[[i]]$R)
    if (length(bad)) {
      tmp <- X[[i]]$L[bad]; X[[i]]$L[bad] <- X[[i]]$R[bad]; X[[i]]$R[bad] <- tmp
    }
  }
  N <- sum(n_i)

  # ---- Precompute Y_ij and group summaries from ORIGINAL data ----
  # Y_ij = (m_ij - mbar_i)^2 + theta * (s_ij - sbar_i)^2
  Ybar_i <- numeric(k)   # \overline Y_{i·}
  den_i  <- numeric(k)   # (1/n_i) \sum_j (Y_ij - \overline Y_{i·})^2
  mu_m   <- mu_s <- numeric(k)  # \bar{X}^C_{i·}, \bar{X}^R_{i·}

  for (i in seq_len(k)) {
    m <- (X[[i]]$L + X[[i]]$R) / 2
    s <- (X[[i]]$R - X[[i]]$L) / 2
    mbar <- mean(m); sbar <- mean(s)
    mu_m[i] <- mbar; mu_s[i] <- sbar
    Y <- (m - mbar)^2 + theta * (s - sbar)^2
    Ybar_i[i] <- mean(Y)
    den_i[i]  <- mean( (Y - mean(Y))^2 )  # equals (1/n_i) * sum_j (...)
  }
  Ydotdot <- sum(n_i * Ybar_i) / N       # \overline Y_{··}

  # ---- Eq. (T_obs): ALWAYS the same observed statistic, independent of `version` ----
  # T_obs =  [ sum_i n_i (Ybar_i - Ydotdot)^2 ] / [ sum_i (1/n_i) sum_j (Y_ij - Ybar_i)^2 ]
  T_obs_num <- sum( n_i * (Ybar_i - Ydotdot)^2 )
  T_obs_den <- sum( den_i )
  T_obs     <- T_obs_num / T_obs_den

  # For Eq. (T1*), the denominator is the ORIGINAL denominator (same as T_obs_den)
  T1star_den <- T_obs_den

  # ---- Bootstrap loop: Steps 2–3 ----
  T_boot <- numeric(B)
  for (b in seq_len(B)) {

    # Step 2: resample intervals within each group and form Y*_ij centered at ORIGINAL (mu_m, mu_s)
    Ybar_star_i <- numeric(k)  # \overline Y^{*}_{i·}
    den_star_T2 <- 0           # denominator for Eq. (T2*)

    for (i in seq_len(k)) {
      idx <- sample.int(n_i[i], n_i[i], replace = TRUE)
      Lst <- X[[i]]$L[idx]; Rst <- X[[i]]$R[idx]
      swap <- which(Lst > Rst)
      if (length(swap)) { tmp <- Lst[swap]; Lst[swap] <- Rst[swap]; Rst[swap] <- tmp }

      m_st <- (Lst + Rst) / 2
      s_st <- (Rst - Lst) / 2

      # Y*_{ij} = (m*_{ij} - mu_m[i])^2 + theta (s*_{ij} - mu_s[i])^2
      Yst <- (m_st - mu_m[i])^2 + theta * (s_st - mu_s[i])^2
      Ybar_star_i[i] <- mean(Yst)

      if (identical(version, "T2")) {
        # Add the i-th contribution to the T2* denominator:
        # (1 / (n_i * Ybar_i^2)) * sum_j (Y*_{ij} - Ybar*_{i·})^2
        den_star_T2 <- den_star_T2 +
          (1 / (n_i[i] * pmax(Ybar_i[i]^2, .Machine$double.eps))) *
          sum( (Yst - mean(Yst))^2 )
      }
    }

    # Step 3: bootstrap statistic numerator (common structure), with version-specific weights
    # delta_i = Ybar*_{i·} - Ybar_{i·}
    delta_i <- Ybar_star_i - Ybar_i

    if (identical(version, "T1")) {
      # Eq. (T1*): w_i = Ydotdot / Ybar_i
      w_i   <- Ydotdot / pmax(Ybar_i, .Machine$double.eps)
      grand <- (1 / N) * sum( n_i * w_i * delta_i )
      num   <- sum( n_i * (w_i * delta_i - grand)^2 )
      denom <- T1star_den
    } else {
      # Eq. (T2*): w_i = 1 / Ybar_i
      w_i   <- 1 / pmax(Ybar_i, .Machine$double.eps)
      grand <- (1 / N) * sum( n_i * w_i * delta_i )
      num   <- sum( n_i * (w_i * delta_i - grand)^2 )
      denom <- den_star_T2
    }

    T_boot[b] <- num / denom
  }

  # ---- Return: T_obs, bootstrap draws, and p-value ----
  list(stat = T_obs, p.value = mean(T_boot >= T_obs), draws = T_boot)
}

