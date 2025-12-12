# ============================================================
# Interval–Levene: Publication Workflow (power/size/timing)
# Compare paper bootstrap-T vs. limiting-F (Alg. 1)
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(future); library(future.apply)
  library(ggplot2); library(tidyr); library(data.table); library(glue)
  library(intervalLevene)
  library(grid)  # <- for unit()
})


# -----------------------------
# 2) RNG + parallel
# -----------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(20250917)
plan(multisession, workers = max(1L, parallel::detectCores() - 2L))


# -----------------------------
# 3) Experiment settings
# -----------------------------
DIST <- "poisson"   # "gamma" or "poisson"
alpha        <- 0.05
SMOKE        <- TRUE
R_boot       <- if (SMOKE)  2000 else 2000
N_reps       <- if (SMOKE)  100 else 10000
tau_grid     <- c(1, 1.2,1.4,1.6,1.8, 2)
target_grp   <- "last"
effect       <- "radius"
rho          <- 0.5
omega        <- 1
tab_dir <- "tables_pub"
fig_dir <- "figs_sim"
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 4) Designs (5 groups, approx-balanced)
# -----------------------------
mk_vec <- function(n) c(n, n+1, n, n+1, n)
approx_vals <- c(10, 30, 50, 70)
designs <- setNames(lapply(approx_vals, mk_vec), paste("Approximately", approx_vals))
stopifnot(anyDuplicated(names(designs)) == 0)


design_levels <- names(designs)
cols <- c("bootstrap T1" = "#56B4E9",
          "limiting (F)" = "#F0E442",
          "bootstrap T2" = "#0072B2")


# -----------------------------
# 5) Helpers
# -----------------------------
.as_df_with_group <- function(sim_list) {
  df <- data.table::rbindlist(
    lapply(seq_along(sim_list), function(i) transform(sim_list[[i]], group = i)),
    use.names = TRUE, fill = TRUE
  )
  df$group <- factor(df$group)
  df
}

# Robust converter: accepts (L,R) OR (C,R) / (center,radius)
.df_to_Xlist <- function(df) {
  stopifnot(is.data.frame(df))
  if (!"group" %in% names(df)) stop("df must contain a 'group' column")

  if (all(c("L","R") %in% names(df))) {
    LR <- df[, c("L","R","group")]
  } else if (all(c("center","radius") %in% names(df))) {
    LR <- transform(
      df,
      L = as.numeric(center) - as.numeric(radius),
      R = as.numeric(center) + as.numeric(radius)
    )[, c("L","R","group")]
  } else {
    stop("df must have either (L,R,group) or (center,radius,group).")
  }

  # (optional but nice) ensure group is factor and fix any L>R
  LR$group <- factor(LR$group)
  bad <- with(LR, L > R)
  if (any(bad)) {
    tmp <- LR$L[bad]; LR$L[bad] <- LR$R[bad]; LR$R[bad] <- tmp
  }

  split(LR[, c("L","R")], LR$group)
}

.sim_once_hetero <- function(gs, rho, tau,
                             dist = c("gamma","poisson"),
                             meanC = 5, scaleC = 2, meanR = 10, scaleR = 1,
                             target = c("last","first","index"),
                             index = NULL,
                             effect = c("radius","center","both")) {
  dist   <- match.arg(dist)
  target <- match.arg(target)
  effect <- match.arg(effect)

  k   <- length(gs)
  idx <- switch(target, last = k, first = 1, index = { if (is.null(index)) k else as.integer(index) })

  lapply(seq_len(k), function(i) {
    if (dist == "gamma") {
      scC <- scaleC; scR <- scaleR
      mnC <- meanC;  mnR <- meanR
      if (i == idx) {
        if (effect %in% c("center","both")) scC <- scaleC * tau
        if (effect %in% c("radius","both")) scR <- scaleR * tau
      }
      out <- try(
        intervalLevene::sim_data(
          group_sizes = gs[i], dist = "gamma", rho = rho,
          center_par = list(mean = mnC, scale = scC),
          radius_par = list(mean = mnR, scale = scR)
        ),
        silent = TRUE
      )
      if (inherits(out, "try-error")) {
        out <- intervalLevene::sim_data(
          group_sizes = gs[i], dist = "gamma", corr = rho,
          params1 = list(mean = mnC, scale = scC),
          params2 = list(mean = mnR, scale = scR)
        )
      }
    } else {
      # --- Poisson: τ multiplies means (lambdas) ---
      lamC <- meanC; lamR <- meanR
      if (i == idx) {
        if (effect %in% c("center","both")) lamC <- meanC * tau
        if (effect %in% c("radius","both")) lamR <- meanR * tau
      }
      out <- try(
        intervalLevene::sim_data(
          group_sizes = gs[i], dist = "poisson", rho = rho,
          center_par = list(mean = lamC),
          radius_par = list(mean = lamR)
        ),
        silent = TRUE
      )
      if (inherits(out, "try-error")) {
        out <- intervalLevene::sim_data(
          group_sizes = gs[i], dist = "poisson", corr = rho,
          params1 = list(mean = lamC),
          params2 = list(mean = lamR)
        )
      }
    }

    d <- out[[1]]
    # Ensure we have both (center,radius) and (L,R)
    nml <- tolower(names(d))
    if (!all(c("center","radius") %in% nml)) {
      cidx <- match(TRUE, nml %in% c("c","center","centre","mid","mean"))
      ridx <- match(TRUE, nml %in% c("r","radius","rad"))
      stopifnot(!is.na(cidx), !is.na(ridx))
      d$center <- as.numeric(d[[cidx]])
      d$radius <- pmax(as.numeric(d[[ridx]]), .Machine$double.eps)
    }
    d$L <- d$center - d$radius
    d$R <- d$center + d$radius
    d$group <- i
    d[, c("center","radius","L","R","group")]
  })
}


crit_lim <- function(k, alpha = 0.05) qchisq(1 - alpha, df = k - 1) / (k - 1)
# -----------------------------
# 6) Single replicate: limiting F, bootstrap T1*, T2*  (granular timing)
# -----------------------------
one_rep_all <- function(gs, rho, tau, alpha = 0.05, R_boot = 1000,
                        target_grp = "last", effect = "radius", omega = 1) {

  # ---- Warm-up to reduce first-call overhead (very short) ----
  # (Only triggers for the first few replicas per worker; negligible cost)
  if (!exists(".il_warmed_up", inherits = FALSE)) {
    invisible(try(stats::runif(10), silent = TRUE))
    .il_warmed_up <<- TRUE
  }

  target_mode <- if (is.numeric(target_grp)) "index" else target_grp
  target_idx  <- if (is.numeric(target_grp)) as.integer(target_grp) else NULL
  k <- length(gs)

  # --- 1) Simulate once ---
  t_sim <- system.time({
    sim_list <- .sim_once_hetero(gs, rho, tau,
                                 dist   = DIST,
                                 target = target_mode, index = target_idx,
                                 effect = effect)
  })[["elapsed"]]

  # --- 2) Assemble df + prep X_list (neutral to all methods) ---
  t_prep <- system.time({
    df     <- .as_df_with_group(sim_list)
    X_list <- .df_to_Xlist(df)
  })[["elapsed"]]

  # --- 3) Compute F-statistic ---
  t_Fstat <- system.time({
    F_obs <- intervalLevene::levene_interval_F(df, omega = omega)
  })[["elapsed"]]

  # --- 4) Limiting decision: just threshold compare ---
  t_F_compare <- system.time({
    q_lim   <- crit_lim(k, alpha)
    rej_lim <- as.integer(F_obs > q_lim)
  })[["elapsed"]]

  # --- 5) Bootstrap T1: draws vs quantile+compare split ---
  t_T1_draws <- NA_real_; t_T1_quant <- NA_real_; rej_T1 <- NA_real_
  tmp <- system.time({
    resT1  <- intervalLevene::interval_bootstrap_RL2012(
      X = X_list, theta = omega, B = R_boot, version = "T1", seed = NULL
    )
  })
  t_T1_draws <- tmp[["elapsed"]]
  t_T1_quant <- system.time({
    q_T1   <- as.numeric(quantile(resT1$draws, 1 - alpha, names = FALSE))
    rej_T1 <- as.integer(resT1$stat > q_T1)
  })[["elapsed"]]

  # --- 6) Bootstrap T2: draws vs quantile+compare split ---
  t_T2_draws <- NA_real_; t_T2_quant <- NA_real_; rej_T2 <- NA_real_
  tmp2 <- system.time({
    resT2  <- intervalLevene::interval_bootstrap_RL2012(
      X = X_list, theta = omega, B = R_boot, version = "T2", seed = NULL
    )
  })
  t_T2_draws <- tmp2[["elapsed"]]
  t_T2_quant <- system.time({
    q_T2   <- as.numeric(quantile(resT2$draws, 1 - alpha, names = FALSE))
    rej_T2 <- as.integer(resT2$stat > q_T2)
  })[["elapsed"]]

  # --- Totals (apples-to-apples) ---
  t_lim_total     <- t_sim + t_prep + t_Fstat + t_F_compare
  t_boot_T1_total <- t_sim + t_prep + t_Fstat + t_T1_draws + t_T1_quant
  t_boot_T2_total <- t_sim + t_prep + t_Fstat + t_T2_draws + t_T2_quant

  setNames(
    c(
      rej_lim          = rej_lim,
      rej_boot_T1      = rej_T1,
      rej_boot_T2      = rej_T2,

      # detailed timings
      t_sim            = t_sim,
      t_prep           = t_prep,
      t_Fstat          = t_Fstat,
      t_F_compare      = t_F_compare,
      t_T1_draws       = t_T1_draws,
      t_T1_quant       = t_T1_quant,
      t_T2_draws       = t_T2_draws,
      t_T2_quant       = t_T2_quant,

      # totals (for fairness)
      t_lim_total      = t_lim_total,
      t_boot_T1_total  = t_boot_T1_total,
      t_boot_T2_total  = t_boot_T2_total,

      F_obs            = F_obs
    ),
    c("rej_lim","rej_boot_T1","rej_boot_T2",
      "t_sim","t_prep","t_Fstat","t_F_compare",
      "t_T1_draws","t_T1_quant","t_T2_draws","t_T2_quant",
      "t_lim_total","t_boot_T1_total","t_boot_T2_total",
      "F_obs")
  )
}

# Always return the same named vector to keep parallel apply stable
safe_one_rep <- function(...) {
  out <- try(one_rep_all(...), silent = TRUE)
  if (inherits(out, "try-error")) {
    setNames(rep(NA_real_, 15),
             c("rej_lim","rej_boot_T1","rej_boot_T2",
               "t_sim","t_prep","t_Fstat","t_F_compare",
               "t_T1_draws","t_T1_quant","t_T2_draws","t_T2_quant",
               "t_lim_total","t_boot_T1_total","t_boot_T2_total","F_obs"))
  } else out
}

# -----------------------------
# 7) Run a full design across tau_grid (F vs T1 vs T2)
# -----------------------------
power_for_design <- function(gs, design_name, tau_grid, N_reps, rho, alpha, R_boot,
                             target_grp, effect, omega) {
  do.call(rbind, lapply(tau_grid, function(tau) {
    reps <- future_sapply(
      X = seq_len(N_reps),
      FUN = function(i, gs, rho, tau, alpha, R_boot, target_grp, effect, omega) {
        safe_one_rep(gs, rho, tau, alpha, R_boot, target_grp, effect, omega)
      },
      gs = gs, rho = rho, tau = tau, alpha = alpha, R_boot = R_boot,
      target_grp = target_grp, effect = effect, omega = omega,
      future.seed = TRUE, simplify = "array",
      future.packages = c("intervalLevene"),
      future.globals = c(  "DIST",
                           "one_rep_all","safe_one_rep",".df_to_Xlist",
                           ".as_df_with_group",".sim_once_hetero","crit_lim"
      )
    )

    m_rej_lim   <- mean(reps["rej_lim",           ], na.rm = TRUE)
    m_rej_T1    <- mean(reps["rej_boot_T1",       ], na.rm = TRUE)
    m_rej_T2    <- mean(reps["rej_boot_T2",       ], na.rm = TRUE)

    # --- choose which timing you want to report into power_tab ---
    # Option A: TOTAL decision time = F_stat + method work  (recommended, apples-to-apples)
    m_t_lim     <- mean(reps["t_lim_total",       ], na.rm = TRUE)
    m_t_T1      <- mean(reps["t_boot_T1_total",   ], na.rm = TRUE)
    m_t_T2      <- mean(reps["t_boot_T2_total",   ], na.rm = TRUE)

    rbind(
      data.frame(
        design      = design_name, k = length(gs), tau = tau,
        method      = "limiting (F)",
        reject_rate = m_rej_lim,
        mean_time_s = m_t_lim,
        N_reps = N_reps, R_boot = R_boot
      ),
      data.frame(
        design      = design_name, k = length(gs), tau = tau,
        method      = "bootstrap T1",
        reject_rate = m_rej_T1,
        mean_time_s = m_t_T1,
        N_reps = N_reps, R_boot = R_boot
      ),
      data.frame(
        design      = design_name, k = length(gs), tau = tau,
        method      = "bootstrap T2",
        reject_rate = m_rej_T2,
        mean_time_s = m_t_T2,
        N_reps = N_reps, R_boot = R_boot
      )
    )
  }))
}


# -----------------------------
# 8) RUN: full grid (power_tab)
# -----------------------------
t_all <- system.time({
  power_tab <- do.call(rbind, lapply(names(designs), function(nm) {
    power_for_design(designs[[nm]], nm, tau_grid, N_reps, rho, alpha, R_boot,
                     target_grp, effect, omega)
  }))
})[["elapsed"]]

# --- after power_tab is built ---
power_tab$method <- factor(power_tab$method, levels = names(cols))





# -----------------------------
# 9) SUMMARIES: size, power, timing (+ CIs)
# -----------------------------
size_tab    <- dplyr::filter(power_tab, tau == 1)
power_only  <- dplyr::filter(power_tab, tau > 1)

# CI helper (per-row N_reps, robust to edge cases)
add_ci <- function(df, nreps_col = "N_reps", rate_col = "reject_rate") {
  df |>
    dplyr::mutate(
      se = sqrt(pmax(0, .data[[rate_col]] * (1 - .data[[rate_col]]) / .data[[nreps_col]])),
      lo = pmax(0, .data[[rate_col]] - 1.96 * se),
      hi = pmin(1, .data[[rate_col]] + 1.96 * se)
    )
}



size_by_design <- size_tab |>
  dplyr::group_by(design, method) |>
  dplyr::summarise(
    reject_rate = mean(reject_rate, na.rm = TRUE),
    N_reps      = dplyr::first(N_reps),
    .groups     = "drop"
  ) |>
  add_ci() |>
  dplyr::mutate(design = factor(design, levels = design_levels))

size_overall <- size_by_design |>
  dplyr::group_by(method) |>
  dplyr::summarise(
    mean_size = mean(reject_rate, na.rm = TRUE),
    sd_size   = sd(reject_rate,   na.rm = TRUE),
    .groups   = "drop"
  )

power_by_tau <- power_only |>
  dplyr::group_by(method, tau) |>
  dplyr::summarise(
    reject_rate = mean(reject_rate, na.rm = TRUE),
    N_reps      = dplyr::first(N_reps),
    .groups     = "drop"
  ) |>
  add_ci()

timing_overall <- power_tab |>
  dplyr::group_by(method) |>
  dplyr::summarise(
    mean_time   = mean(mean_time_s, na.rm = TRUE),
    median_time = stats::median(mean_time_s, na.rm = TRUE),
    .groups     = "drop"
  )

speedup_tab <- timing_overall |>
  dplyr::filter(method %in% c("bootstrap T1", "bootstrap T2", "limiting (F)")) |>
  tidyr::pivot_wider(
    names_from = method,
    values_from = mean_time,
    names_repair = "minimal",
    values_fill = NA_real_
  )

# Compute ratios using .data[[...]] to safely reference names with spaces
speedup_tab <- speedup_tab |>
  dplyr::mutate(
    speedup_T1_over_lim = .data[["bootstrap T1"]] / .data[["limiting (F)"]],
    speedup_T2_over_lim = .data[["bootstrap T2"]] / .data[["limiting (F)"]]
  )


# Optional console line:
if (nrow(speedup_tab)) {
  cat(glue("\nSpeedup (mean time): T1/lim={round(speedup_tab$speedup_T1_over_lim,2)}x; ",
           "T2/lim={round(speedup_tab$speedup_T2_over_lim,2)}x\n"))
}

readr::write_csv(size_overall,   file.path(tab_dir, glue("size_overall_R{R_boot}_N{N_reps}.csv")))
readr::write_csv(power_by_tau,   file.path(tab_dir, glue("power_by_tau_R{R_boot}_N{N_reps}.csv")))
readr::write_csv(timing_overall, file.path(tab_dir, glue("timing_overall_R{R_boot}_N{N_reps}.csv")))
# -----------------------------
# 10) FIGURES (saved) with CIs
# -----------------------------
# Type I error bars by design


size_by_design$method <- factor(size_by_design$method, levels = names(cols))

p_size <- ggplot(size_by_design, aes(x = design, y = reject_rate, fill = method)) +
  geom_hline(yintercept = alpha, linetype = 2) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.95,
           show.legend = TRUE) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.7),
                width = 0.2, linewidth = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = function(x) gsub("^\\s*Approximately\\s+", "", x)) +
  labs(y = "Empirical size", x = "Groups zizes") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "top",         # <- put legend at bottom
    legend.title = element_blank()
  )
p_size

ggsave(file.path(fig_dir, glue::glue("type1erro_vs_tau_by_design_{DIST}_R{R_boot}_N{N_reps}.png")),
       p_size, width = 9, height = 6.5, dpi = 600)

# --- Power vs tau, FACETTED by design (with CIs) ------------------------------
power_by_tau_design <- power_only |>
  dplyr::group_by(design, method, tau) |>
  dplyr::summarise(
    reject_rate = mean(reject_rate, na.rm = TRUE),
    N_reps      = dplyr::first(N_reps),
    .groups     = "drop"
  ) |>
  add_ci() |>
  dplyr::mutate(
    design = factor(design, levels = design_levels),
    method = factor(method, levels = names(cols))
  )

p_power_tau_facet <- ggplot(power_by_tau_design,
                            aes(x = tau, y = reject_rate, color = method, fill = method)) +
  # geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.9, stroke = 0.4) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = sort(unique(power_by_tau_design$tau))) +
  labs(
    x = expression(paste("Heterogeneity factor  ", tau)),
    y = "Power (mean rejection rate)"
  ) +
  facet_wrap(~ design, ncol = 2) +
  theme_minimal() +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(0.5, "lines"),
    strip.text       = element_text(size = 16, face = "bold", lineheight = 0.9,
                                    margin = margin(t = 2, r = 6, b = 2, l = 6)),
    strip.background = element_rect(fill = "gray96", color = "gray85", linewidth = 0.6),
    axis.title       = element_text(size = 16, face = "bold"),
    axis.text        = element_text(size = 16),
    panel.border     = element_rect(color = "gray75", fill = NA, linewidth = 1.4),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.8),
    panel.grid.minor = element_line(color = "gray96", linewidth = 0.8),
    plot.margin      = margin(15, 15, 15, 15),
    legend.position  = "top",
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 14)
  )


p_power_tau_facet
ggsave(file.path(fig_dir, glue::glue("power_vs_tau_by_design_{DIST}_R{R_boot}_N{N_reps}.png")),
       p_power_tau_facet, width = 9, height = 6.5, dpi = 600)
#===============================================================================
# Build timing summary per design × tau × method
timing_by_tau_design <- power_tab |>
  dplyr::group_by(design, method, tau) |>
  dplyr::summarise(
    mean_time_s = mean(mean_time_s, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    design = factor(design, levels = design_levels),
    method = factor(method, levels = names(cols))
  )

# Make sure the factor level order matches the palette
timing_by_tau_design$method <- factor(timing_by_tau_design$method, levels = names(cols))
timing_by_tau_design$design <- factor(timing_by_tau_design$design, levels = names(designs))




# 1) Summarize over tau -> one point per (design, method)
df_sum <- timing_by_tau_design %>%
  mutate(
    design_num = dplyr::coalesce(
      suppressWarnings(as.numeric(design)),
      readr::parse_number(as.character(design))
    )
  ) %>%
  group_by(method, design, design_num) %>%
  summarise(
    mean_time = mean(mean_time_s, na.rm = TRUE),
    sd_time   = sd(mean_time_s, na.rm = TRUE),
    n_tau     = dplyr::n(),
    se_time   = sd_time / sqrt(pmax(n_tau, 1)),
    .groups   = "drop"
  ) %>%
  arrange(design_num)

# 2) Nice labels for x-axis if your design column had "Approximately XX"
x_breaks <- sort(unique(df_sum$design_num))
x_labels <- if (any(is.na(x_breaks))) unique(df_sum$design) else x_breaks

# ensure same legend order + colors as other plots
df_sum$method <- factor(df_sum$method, levels = names(cols))

# 3) Plot
p_time_by_design <- ggplot(
  df_sum,
  aes(x = design_num, y = mean_time, color = method, group = method)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(mean_time - se_time, 0),
                    ymax = mean_time + se_time),
                width = 0.8, alpha = 0.5) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_color_manual(values = cols) +        # ▶ use same palette
  labs(

    x = "Design number",
    y = "Mean runtime per replicate (seconds)",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    axis.title       = element_text(size = 16, face = "bold"),
    axis.text        = element_text(size = 14),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.8),
    panel.grid.minor = element_line(color = "gray96", linewidth = 0.6),
    panel.border     = element_rect(color = "gray75", fill = NA, linewidth = 1.2),
    legend.position  = "top",
    plot.margin      = margin(15, 15, 15, 15)
  )


p_time_by_design

# 4) Save
ggsave(
  file.path(fig_dir, glue::glue("runtime_by_design_{DIST}_R{R_boot}_N{N_reps}.png")),
  p_time_by_design, width = 8.5, height = 5.5, dpi = 600
)

# -----------------------------
# 11) TABLE 1: q_emp, q_lim, c_alpha, exceedance (τ = 1) based on F
# -----------------------------
null_F_draws <- do.call(rbind, lapply(names(designs), function(nm) {
  reps <- future_sapply(
    X = seq_len(N_reps),
    FUN = function(i, gs, rho, alpha, R_boot, target_grp, effect, omega) {
      out <- try(
        one_rep_all(gs, rho, tau = 1, alpha = alpha, R_boot = R_boot,
                    target_grp = target_grp, effect = effect, omega = omega),
        silent = TRUE
      )
      if (inherits(out, "try-error")) NA_real_ else out["F_obs"]
    },
    gs = designs[[nm]], rho = rho, alpha = alpha, R_boot = R_boot,
    target_grp = target_grp, effect = effect, omega = omega,
    future.seed = TRUE, simplify = "array",
    future.packages = c("intervalLevene"),
    future.globals = c( "DIST",
                        "one_rep_all",".df_to_Xlist",
                        ".as_df_with_group",".sim_once_hetero","crit_lim"
    )
  )
  data.frame(design = nm, F_obs = as.numeric(reps))
}))


# Save raw null draws
write_csv(null_F_draws, file.path(tab_dir, glue("null_F_draws_R{R_boot}_N{N_reps}.csv")))

# All designs use k = length(gs); if yours vary by design, compute per design.
k_fixed <- length(designs[[1]])
q_lim_0_95 <- qchisq(1 - alpha, df = k_fixed - 1) / (k_fixed - 1)

tab1 <- null_F_draws %>%
  group_by(design) %>%
  summarise(
    q_emp_0_95      = as.numeric(quantile(F_obs, 0.95, na.rm = TRUE, names = FALSE)),
    q_lim_0_95      = q_lim_0_95,
    c_alpha         = q_emp_0_95 / q_lim_0_95,
    exceed_at_q_lim = mean(F_obs > q_lim_0_95, na.rm = TRUE),
    n_draws         = sum(!is.na(F_obs)),
    .groups = "drop"
  ) %>%
  arrange(factor(design, levels = names(designs)))

# Save CSV and TeX preview
write_csv(tab1, file.path(tab_dir, glue("table1_critvals_R{R_boot}_N{N_reps}.csv")))

tex_lines <- c(
  "\\begin{tabular}{lrrrr}",
  "\\toprule",
  "Design & $q^{\\text{emp}}_{0.95}$ & $q^{\\text{lim}}_{0.95}$ & $c_\\alpha$ & $\\Pr\\{F>q^{\\text{lim}}_{0.95}\\}$\\\\",
  "\\midrule",
  paste0(
    sprintf("%s & %.4f & %.4f & %.3f & %.3f\\\\",
            tab1$design, tab1$q_emp_0_95, tab1$q_lim_0_95, tab1$c_alpha, tab1$exceed_at_q_lim),
    collapse = "\n"
  ),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(tex_lines, con = file.path(tab_dir, glue("table1_critvals_R{R_boot}_N{N_reps}.tex")))

# -----------------------------
# 12) Console summary (guarded)
# -----------------------------
cat("\n== TABLE 1 (Critical values & exceedance at τ = 1) ==\n"); print(tab1)

# These objects come from earlier sections; print them if they exist.
if (exists("size_overall"))   { cat("\n== SIZE OVERALL ==\n");   print(size_overall) }
if (exists("power_by_tau"))   { cat("\n== POWER BY TAU ==\n");   print(power_by_tau) }
if (exists("timing_overall")) { cat("\n== TIMING OVERALL ==\n"); print(timing_overall) }


cat(glue("\nArtifacts saved in: {normalizePath(tab_dir)} and {normalizePath(fig_dir)}\n"))

# Clean up plan
plan(sequential)
nbrOfWorkers()
