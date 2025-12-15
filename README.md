
<!-- README.md is generated from README.Rmd. Please edit that file -->

## intervalLevene

<!-- badges: start -->
<!-- badges: end -->

intervalLevene provides Levene-type tests for homogeneity of dispersion
with interval-valued data, using a center–radius representation of
intervals. It includes a fast limiting reference distribution for an
ANOVA-style statistic, bootstrap competitors for comparison, and
simulation utilities used in the accompanying manuscript.

\##Installation

You can install the development version from GitHub with:

``` r
## install.packages("pak")
pak::pak("Mahamoussa321/intervalLevene")
```

## Data: Utah SWE case study

The package ships a dataset used in the Utah SWE case study.

``` r
data("soft_data_utah")
str(soft_data_utah)
#> Classes 'tbl_df', 'tbl' and 'data.frame':    1820 obs. of  5 variables:
#>  $ lon     : num  -114 -113 -113 -113 -112 ...
#>  $ lat     : num  42 42 42 42 42 ...
#>  $ KG_major: chr  "B" "B" "B" "D" ...
#>  $ center  : num  3.33e-04 4.51e-03 1.92e-05 1.61e-04 6.84e-05 ...
#>  $ radius  : num  2.354 1.039 0.131 1.136 0.484 ...
table(soft_data_utah$KG_major)
#> 
#>    B    C    D 
#> 1180  140  500
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(intervalLevene)

data("soft_data_utah", package = "intervalLevene")
#> Warning in data("soft_data_utah", package = "intervalLevene"): data set
#> 'soft_data_utah' not found

# Keep B/C/D only
utah <- subset(soft_data_utah, KG_major %in% c("B", "C", "D"))
utah$KG_major <- factor(utah$KG_major, levels = c("B", "C", "D"))

# Interval–Levene statistic (limiting reference)
omega <- 1
F_obs <- levene_interval_F(
  df    = utah[, c("center", "radius")],
  group = utah$KG_major,
  omega = omega
)

# Limiting reference: (k-1) * F  ~  ChiSquare(df = k-1)
k <- nlevels(utah$KG_major)
df1 <- k - 1
p_lim <- 1 - pchisq(df1 * F_obs, df = df1)

# Decision at alpha = 0.05
alpha <- 0.05
crit_lim <- qchisq(1 - alpha, df = df1) / df1
reject <- (F_obs > crit_lim)

list(
  F_obs    = F_obs,
  p_value  = p_lim,
  crit_0_05 = crit_lim,
  reject_0_05 = reject
)
#> $F_obs
#> [1] 15.53574
#> 
#> $p_value
#> [1] 1.790259e-07
#> 
#> $crit_0_05
#> [1] 2.995732
#> 
#> $reject_0_05
#> [1] TRUE
```

## Reproducibility

Paper outputs (tables/figures) are stored in folders such as
tables_pub/, tables_sim/, and figs_sim. These are kept in the GitHub
repository for reproducibility, but are excluded from the package build
via .Rbuildignore.
