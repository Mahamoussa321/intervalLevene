# ============================================================
# Utah case study (B, C, D): Interval–Levene + publication plots
# Uses intervalLevene::levene_interval_F() with limiting reference
# Data source: intervalLevene::soft_data_utah (packaged dataset)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(intervalLevene)   # your package
})

# ----------------------------
# Output & settings
# ----------------------------
OUT_DIR <- "figs_sim"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
OMEGA <- 1
set.seed(42)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ----------------------------
# Load packaged Utah data
# ----------------------------
dat_utah <- intervalLevene::soft_data_utah

# Safety checks (helps if the dataset format changes later)
needed <- c("lon", "lat", "KG_major", "center", "radius")
miss <- setdiff(needed, names(dat_utah))
if (length(miss) > 0) stop("soft_data_utah is missing columns: ", paste(miss, collapse = ", "))

dat_utah <- dat_utah %>%
  filter(!is.na(lon), !is.na(lat),
         !is.na(KG_major),
         !is.na(center),
         !is.na(radius)) %>%
  mutate(
    KG_major = as.character(KG_major),
    center   = as.numeric(center),
    radius   = pmax(0, as.numeric(radius))
  )

# ----------------------------
# Keep only Köppen major groups B, C, D
# + add human-readable labels used everywhere downstream
# ----------------------------
dat_utah_bcd <- dat_utah %>%
  filter(KG_major %in% c("B","C","D")) %>%
  mutate(
    KG_major = factor(KG_major, levels = c("B","C","D")),
    KG_label = factor(KG_major, levels = c("B","C","D"),
                      labels = c("Arid", "Temperate", "Cold"))
  )

# Counts
n_by_major <- dat_utah_bcd %>%
  count(KG_major) %>%
  tidyr::complete(KG_major = factor(c("B","C","D"), levels = c("B","C","D")),
                  fill = list(n = 0))

N <- nrow(dat_utah_bcd)
k <- nlevels(dat_utah_bcd$KG_major)

message(glue("Utah (B,C,D): N = {N}, counts => ",
             paste0(levels(dat_utah_bcd$KG_major), "=", n_by_major$n, collapse = ", ")))

if (k < 2) stop("Need at least two groups after filtering.")

# ----------------------------
# Interval–Levene test (limiting reference)
# ----------------------------
F_obs <- intervalLevene::levene_interval_F(
  df    = dat_utah_bcd[, c("center","radius")],
  group = dat_utah_bcd$KG_major,
  omega = OMEGA
)

df1      <- k - 1L
crit_lim <- qchisq(1 - ALPHA, df = df1) / df1
p_lim    <- 1 - pchisq(df1 * F_obs, df = df1)
reject   <- F_obs > crit_lim

# ----------------------------
# Console report
# ----------------------------
cat("\n================ Interval–Levene (Utah; Arid/Temperate/Cold) ================\n")
cat(glue("Groups (k): {k}\n"))
cat(glue("Total N    : {N}\n"))
cat(glue("Counts     : B={n_by_major$n[n_by_major$KG_major=='B']}, ",
         "C={n_by_major$n[n_by_major$KG_major=='C']}, ",
         "D={n_by_major$n[n_by_major$KG_major=='D']}\n"))
cat(glue("omega      : {OMEGA}\n"))
cat(glue("Statistic  : F = {sprintf('%.6f', F_obs)}\n"))
cat(glue("Cutoff     : {sprintf('%.4f', crit_lim)} at alpha = {ALPHA}\n"))
cat(glue("p-value    : {format(p_lim, digits = 6)}\n"))
cat(glue("Decision   : {ifelse(reject, 'REJECT H0 (heterogeneity detected)',
                               'Fail to reject H0')}\n"))
cat(  "============================================================================\n\n")


# ----------------------------
# Publication figure: Center vs Radius by major group (B,C,D)
# ----------------------------
long_bcd <- dat_utah_bcd %>%
  select(KG_label, center, radius) %>%
  tidyr::pivot_longer(c(center, radius), names_to = "measure", values_to = "value")

p_utah_BCD <- ggplot(long_bcd, aes(x = measure, y = value, fill = measure)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.6, linewidth = 0.3) +
  facet_wrap(~ KG_label, scales = "free_y") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = scales::pretty_breaks(5)
  ) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "", y = "Value") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "utah_BCD.png"), p_utah_BCD, width = 8, height = 5, dpi = 300)

# ==============================================================================
# === Utah map: per-location mean "center", colored by group (Arid/Temperate/Cold)
# ==============================================================================
drop_zero <- FALSE

# Utah boundary
usa_states <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf") |>
  sf::st_transform(4326)
utah_sf <- dplyr::filter(usa_states, name == "Utah")

# Per-location mean center (since we no longer have swe_mean; center is the analogue)
utah_station_means <- dat_utah_bcd %>%
  { if (drop_zero) dplyr::filter(., center > 0) else . } %>%
  mutate(
    Group = factor(KG_major, levels = c("B","C","D"),
                   labels = c("Arid","Temperate","Cold"))
  ) %>%
  group_by(lon, lat, Group) %>%
  summarise(mean_center = mean(center, na.rm = TRUE), .groups = "drop")

utah_station_sf <- sf::st_as_sf(utah_station_means, coords = c("lon","lat"), crs = 4326)

p_utah_means <- ggplot() +
  geom_sf(data = utah_sf, fill = NA, color = "grey20", linewidth = 0.6) +
  geom_sf(data = utah_station_sf,
          aes(color = Group, size = mean_center),
          alpha = 0.75) +
  scale_color_viridis_d(option = "plasma", name = "Major group") +
  scale_size_continuous(name = "Mean center") +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave(file.path(OUT_DIR, "utah_station_mean_center_map.png"),
       p_utah_means, width = 7.5, height = 5.5, dpi = 300)
