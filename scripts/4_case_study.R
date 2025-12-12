# ============================================================
# Utah case study (B, C, D): Interval–Levene + publication plots
# Uses intervalLevene::levene_interval_F() with limiting reference
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(glue)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(intervalLevene)      # your package
})

# ----------------------------
# Paths & settings
# ----------------------------
DATA_FILE <- "data/soft_data_USA.csv"
OUT_DIR   <- "figs_sim"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
OMEGA <- 1

set.seed(42)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ----------------------------
# Load data
# ----------------------------
message("Reading data: ", DATA_FILE)
usa_raw <- read_csv(DATA_FILE, show_col_types = FALSE)

# Build minimal columns and drop rows with missing essentials
dat_us <- usa_raw %>%
  filter(!is.na(lon), !is.na(lat),
         !is.na(KG_major),
         !is.na(swe_mean),
         !is.na(lower_95), !is.na(upper_95)) %>%
  transmute(
    lon, lat,
    KG_major = as.character(KG_major),
    center   = swe_mean,
    radius   = pmax(0, (upper_95 - lower_95)/2)   # ensure nonnegative
  )

# ----------------------------
# Spatial filter: Utah only
# ----------------------------
usa_states <- ne_states(country = "United States of America", returnclass = "sf") |>
  st_transform(4326)
utah_sf <- dplyr::filter(usa_states, name == "Utah")

pts_sf <- st_as_sf(dat_us, coords = c("lon","lat"), crs = 4326, remove = FALSE)
utah_pts <- st_join(pts_sf, utah_sf["name"], join = st_within, left = FALSE) |>
  st_drop_geometry()

# ----------------------------
# Keep only Köppen major groups B, C, D
# + add human-readable labels used everywhere downstream
# ----------------------------
dat_utah_bcd <- utah_pts %>%
  filter(KG_major %in% c("B","C","D")) %>%
  mutate(
    KG_major = factor(KG_major, levels = c("B","C","D")),
    KG_label = factor(KG_major, levels = c("B","C","D"),
                      labels = c("Arid", "Temperate", "Cold"))
  )

# Quick counts (by code and by label)
n_by_major <- dat_utah_bcd %>% count(KG_major) %>%
  tidyr::complete(KG_major = factor(c("B","C","D"), levels = c("B","C","D")),
                  fill = list(n = 0))
n_by_label <- dat_utah_bcd %>% count(KG_label)

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
  group = dat_utah_bcd$KG_major,   # test is by major-code groups (B,C,D)
  omega = OMEGA
)

df1 <- k - 1L
crit_lim <- qchisq(1 - ALPHA, df = df1) / df1         # χ^2_{k-1}/(k-1)
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
# Save artifacts (TXT + CSV)
# ----------------------------
out_txt <- file.path(OUT_DIR, glue("utah_BCD_interval_levene_limit_{stamp}.txt"))
out_csv_sizes <- file.path(OUT_DIR, glue("utah_BCD_group_sizes_{stamp}.csv"))
out_csv_summary <- file.path(OUT_DIR, glue("utah_BCD_test_summary_{stamp}.csv"))

writeLines(
  c(
    "Utah Interval–Levene case study (limiting reference)",
    glue("timestamp : {stamp}"),
    glue("data_file : {DATA_FILE}"),
    glue("k         : {k}"),
    glue("N         : {N}"),
    glue("counts    : B={n_by_major$n[n_by_major$KG_major=='B']}, ",
         "C={n_by_major$n[n_by_major$KG_major=='C']}, ",
         "D={n_by_major$n[n_by_major$KG_major=='D']}"),
    glue("omega     : {OMEGA}"),
    glue("F         : {sprintf('%.6f', F_obs)}"),
    glue("cutoff    : {sprintf('%.4f', crit_lim)} (alpha = {ALPHA})"),
    glue("p_value   : {format(p_lim, digits = 6)}"),
    glue("decision  : {ifelse(reject, 'REJECT H0', 'Fail to reject H0')}")
  ),
  con = out_txt
)

write_csv(n_by_major, out_csv_sizes)

write_csv(
  tibble(
    alpha = ALPHA, omega = OMEGA, k = k, N = N,
    F_obs = F_obs, crit_lim = crit_lim, p_value = p_lim,
    decision = ifelse(reject, "REJECT H0", "Fail to reject H0")
  ),
  out_csv_summary
)

message("Saved:")
message("  - ", out_txt)
message("  - ", out_csv_sizes)
message("  - ", out_csv_summary)

# ----------------------------
# Publication figure: Center vs Radius by major group (B,C,D)
#   - Zeros INCLUDED; pseudo-log y (handles zeros gracefully)
#   - Facet titles use Arid/Temperate/Cold
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
# === Utah map: per-station mean SWE, colored by group (Arid/Temperate/Cold) ===
# ==============================================================================

# Toggle this if you want to exclude zeros before averaging
drop_zero <- FALSE

# Per-station means within Utah, labeled by human-readable group
utah_station_means <- usa_raw %>%
  dplyr::filter(!is.na(lon), !is.na(lat), !is.na(KG_major), !is.na(swe_mean)) %>%
  sf::st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE) %>%
  sf::st_join(utah_sf["name"], join = sf::st_within, left = FALSE) %>%
  sf::st_drop_geometry() %>%
  dplyr::filter(KG_major %in% c("B","C","D")) %>%
  { if (drop_zero) dplyr::filter(., swe_mean > 0) else . } %>%
  dplyr::mutate(
    Group = factor(KG_major, levels = c("B","C","D"),
                   labels = c("Arid","Temperate","Cold"))
  ) %>%
  dplyr::group_by(lon, lat, Group) %>%
  dplyr::summarise(mean_swe = mean(swe_mean, na.rm = TRUE), .groups = "drop")

utah_station_sf <- sf::st_as_sf(utah_station_means, coords = c("lon","lat"), crs = 4326)

p_utah_means <- ggplot() +
  geom_sf(data = utah_sf, fill = NA, color = "grey20", linewidth = 0.6) +
  geom_sf(data = utah_station_sf,
          aes(color = Group, size = mean_swe),
          alpha = 0.75) +
  scale_color_viridis_d(option = "plasma", name = "Major group") +
  scale_size_continuous(name = "Mean SWE (mm)") +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave(file.path(OUT_DIR, "utah_station_mean_SWE_map.png"),
       p_utah_means, width = 7.5, height = 5.5, dpi = 300)
