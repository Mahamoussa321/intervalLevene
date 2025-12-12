
#===============================================================================
library(intervalLevene)
library(future.apply); plan(multisession, workers = max(1, parallel::detectCores()-1))
library(dplyr); library(ggplot2)

# Designs you want to validate
designs <- list(

  "Approximately 10" = c(10,11,10,11,10),
  #"Approximately 20" = c(20,21,20,21,20),
  "Approximately 30" = c(30,31,30,31,30),
  #"Approximately 40" = c(40,41,40,41,40),
  "Approximately 50" = c(50,51,50,51,50),
  #"Approximately 60" = c(60,61,60,61,60))
  "Approximately 70" = c(70,71,70,71,70))
# "Approximately 80" = c(80,81,80,81,80),
# "Approximately 90" = c(90,91,90,91,90))

# "Approximately 100" = c(100,101,100,101,100),
# "Approximately 110" = c(110,111,110,111,110),
# "Approximately 120" = c(120,121,120,121,120))




B <- 10000   # empirical draws (increase if you like)
M <- 10000  # limit draws
rho <- 0.5
omega <- 1

RNGkind("L'Ecuyer-CMRG")  # good for parallel + reproducible futures
set.seed(2555)
# helper: one empirical F draw for a given design
one_draw <- function(gs){
  sim <- intervalLevene::sim_data(gs, dist = "poisson",
                                  center_par = list(mean=5, disp=2),
                                  radius_par = list(mean=10, disp=1),
                                  rho = rho)
  df <- data.table::rbindlist(sim); df$group <- factor(df$group)
  intervalLevene::levene_interval_F(df, df$group, omega = omega)
}

# run everything (parallel)
emp <- lapply(designs, function(gs) {
  future_sapply(seq_len(B), function(i) one_draw(gs), future.seed = TRUE)
})
lim <- lapply(designs, function(gs) stats::rchisq(M, df = length(gs)-1)/(length(gs)-1))

# plot (faceted overlay)
# --- build frames with common labels ---

####################################################################################

# Desired facet order
level_order <- c( "Approximately 10", "Approximately 30", "Approximately 50",
                  #"Approximately 50")
                  #"Approximately 60")
                  "Approximately 70")
# "Approximately 80", "Approximately 90")
# # "Approximately 100",
# "Approximately 110", "Approximately 120")

# Build frames with panel as an ordered factor
emp_df <- dplyr::bind_rows(lapply(names(emp), function(nm)
  data.frame(panel = factor(nm, levels = level_order), value = emp[[nm]])
))

lim_df <- dplyr::bind_rows(lapply(names(lim), function(nm)
  data.frame(panel = factor(nm, levels = level_order), value = lim[[nm]])
))

q95_emp_df <- emp_df |>
  dplyr::group_by(panel) |>
  dplyr::summarise(q95_emp = as.numeric(stats::quantile(value, 0.95)), .groups = "drop")

q95_lim_df <- data.frame(
  panel   = factor(names(designs), levels = level_order),
  q95_lim = sapply(designs, function(gs) stats::qchisq(0.95, df = length(gs)-1)/(length(gs)-1))
)


## Plot — single legend driven by the vlines only ----
p <- ggplot() +
  # histogram: constant colors, NO legend
  geom_histogram(
    data = emp_df,
    aes(value, y = after_stat(density)),
    bins = 70, color ="#F7A5A5" , fill ="#FFDBB6",
    show.legend = FALSE
  ) +
  # limiting density curve: constant color, NO legend
  geom_density(
    data = lim_df,
    aes(value),
    color = "#5D688A", linewidth = 0.9,
    show.legend = FALSE
  ) +
  # 95th percentiles: these TWO layers create the ONLY legend
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
    values = c(Data = "#FFDBB6", Limit = "#5D688A")
  ) +
  guides(
    # make legend show line samples with the right linetypes
    color = guide_legend(override.aes = list(linetype = c("dashed", "solid"), linewidth = 0.9))
  ) +
  labs(x = "value", y = "density") +
  theme_bw(12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85), color = NA)
  )

print(p)
dir.create("figs_sim", showWarnings = FALSE)
ggsave("figs_sim/interval_levene_facets_poisson.png",
       p, width = 11, height = 8, dpi = 300)

####################################################################################

# numeric 95th percentiles (table for the paper)
qtab <- data.frame(
  panel = names(designs),
  k     = sapply(designs, length),
  q95_emp = sapply(emp, \(x) as.numeric(quantile(x, .95))),
  q95_lim = sapply(designs, \(gs) stats::qchisq(.95, df=length(gs)-1)/(length(gs)-1))
)
readr::write_csv(qtab, "tables_sim/interval_levene_q95.csv")
qtab

