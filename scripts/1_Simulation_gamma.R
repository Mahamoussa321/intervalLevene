# ============================================================
# Interval-Levene validation: empirical vs limiting reference
# ============================================================

suppressPackageStartupMessages({
  library(intervalLevene)
  library(future)         # plan() lives here
  library(future.apply)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# Parallel plan (use all-1 cores)
# -----------------------------
future::plan(multisession, workers = max(1, parallel::detectCores() - 1))

# -----------------------------
# OUTPUT PATHS
# -----------------------------
fig_dir <- "figs_sim"
tab_dir <- "tables_sim"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# DESIGNS (vectors of group sizes)
# NOTE: Labels now match sizes (~10, ~20, ..., ~90).
# -----------------------------
designs <- list(
  "Approximately 10" = c(10, 11, 10, 11, 10),
  "Approximately 20" = c(20, 21, 20, 21, 20),
  "Approximately 30" = c(30, 31, 30, 31, 30),
  "Approximately 40" = c(40, 41, 40, 41, 40),
  "Approximately 50" = c(50, 51, 50, 51, 50),
  "Approximately 60" = c(60, 61, 60, 61, 60),
  "Approximately 70" = c(70, 71, 70, 71, 70),
  "Approximately 80" = c(80, 81, 80, 81, 80),
  "Approximately 90" = c(90, 91, 90, 91, 90)
)

# -----------------------------
# SIMULATION CONTROLS
# -----------------------------
B     <- 10000   # Monte Carlo (empirical) draws PER DESIGN
M     <- 10000   # Limiting-law draws PER DESIGN
rho   <- 0.5     # center–radius copula correlation in sim_data()
omega <- 1       # weight on radius term in Levene Y

# Reproducible parallel RNG
RNGkind("L'Ecuyer-CMRG")
set.seed(258)

# -----------------------------
# One empirical F draw for a design 'gs'
# -----------------------------
one_draw <- function(gs) {
  sim <- intervalLevene::sim_data(
    group_sizes = gs,
    dist        = "gamma",
    center_par  = list(mean = 5,  scale = 2),
    radius_par  = list(mean = 10, scale = 1),
    rho         = rho
  )
  df <- data.table::rbindlist(sim)
  df$group <- factor(df$group)
  intervalLevene::levene_interval_F(df, df$group, omega = omega)
}

# -----------------------------
# Run simulations (parallel)
# emp: empirical distribution via data simulation + Levene F
# lim: theoretical limiting reference χ^2_{k-1}/(k-1) (no data simulation)
# -----------------------------
emp <- lapply(designs, function(gs) {
  future.apply::future_sapply(seq_len(B), function(i) one_draw(gs), future.seed = TRUE)
})
lim <- lapply(designs, function(gs) {
  k <- length(gs)
  stats::rchisq(M, df = k - 1) / (k - 1)
})

# -----------------------------
# Faceting order (choose subset or all)
# If you want all, set: level_order <- names(designs)
# -----------------------------
level_order <- c("Approximately 10", "Approximately 30", "Approximately 50", "Approximately 70")
level_order <- intersect(level_order, names(designs))
if (length(level_order) == 0L) level_order <- names(designs)

# -----------------------------
# Build long frames for plotting
# -----------------------------
emp_df <- dplyr::bind_rows(lapply(names(emp), function(nm) {
  data.frame(panel = factor(nm, levels = level_order), value = emp[[nm]])
}))
lim_df <- dplyr::bind_rows(lapply(names(lim), function(nm) {
  data.frame(panel = factor(nm, levels = level_order), value = lim[[nm]])
}))

# 95th percentiles (empirical vs limiting)
q95_emp_df <- emp_df %>%
  dplyr::filter(!is.na(panel)) %>%
  dplyr::group_by(panel) %>%
  dplyr::summarise(q95_emp = as.numeric(stats::quantile(value, 0.95)), .groups = "drop")

q95_lim_df <- data.frame(
  panel   = factor(names(designs), levels = level_order),
  q95_lim = sapply(designs, function(gs) stats::qchisq(0.95, df = length(gs) - 1) / (length(gs) - 1))
) %>% dplyr::filter(!is.na(panel))

# -----------------------------
# Plot (empirical histogram + limiting density + 95th lines)
# -----------------------------
p <- ggplot() +
  geom_histogram(
    data = emp_df %>% dplyr::filter(!is.na(panel)),
    aes(value, y = after_stat(density)),
    bins = 70, color = "#F9EC7E", fill = "#E3CCB2",
    show.legend = FALSE
  ) +
  geom_density(
    data = lim_df %>% dplyr::filter(!is.na(panel)),
    aes(value),
    color = "#E26274", linewidth = 0.9,
    show.legend = FALSE
  ) +
  geom_vline(
    data = q95_emp_df,
    aes(xintercept = q95_emp, color = "Data"),
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_vline(
    data = q95_lim_df,
    aes(xintercept = q95_lim, color = "Limit"),
    linetype = "solid", linewidth = 0.8
  ) +
  facet_wrap(~ panel, scales = "free_x") +
  scale_color_manual(
    name = "95th Quantile",
    values = c(Data = "#E3CCB2", Limit = "#E26274")
  ) +
  guides(color = guide_legend(
    override.aes = list(linetype = c("dashed", "solid"), linewidth = 0.9)
  )) +
  labs(x = "value", y = "density") +
  theme_bw(12) +
  theme(
    strip.text          = element_text(face = "bold"),
    legend.position     = c(0.98, 0.98),
    legend.justification= c(1, 1),
    legend.background   = element_rect(fill = scales::alpha("white", 0.85), color = NA)
  )

print(p)
ggsave(file.path(fig_dir, "interval_levene_facets.png"), p, width = 11, height = 8, dpi = 300)

# -----------------------------
# Numeric table (q95 empirical vs limiting) for the paper
# -----------------------------
qtab <- data.frame(
  panel   = names(designs),
  k       = sapply(designs, length),
  q95_emp = sapply(emp, function(x) as.numeric(stats::quantile(x, 0.95))),
  q95_lim = sapply(designs, function(gs) stats::qchisq(0.95, df = length(gs) - 1) / (length(gs) - 1))
) %>%
  dplyr::mutate(diff_emp_minus_lim = q95_emp - q95_lim)

readr::write_csv(qtab, file.path(tab_dir, "interval_levene_q95.csv"))
qtab
