#' Utah SWE “soft” interval data (case study)
#'
#' A dataset for the Utah case study used in the package examples/scripts.
#' It contains station coordinates, Köppen–Geiger major climate group, and
#' interval representation of SWE via (center, radius).
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{lon}{Longitude in decimal degrees (WGS84).}
#'   \item{lat}{Latitude in decimal degrees (WGS84).}
#'   \item{KG_major}{Köppen–Geiger major group (e.g., B, C, D).}
#'   \item{center}{Mean SWE value.}
#'   \item{radius}{Half-width of the 95% interval for SWE (nonnegative).}
#' }
#' @keywords datasets
"soft_data_utah"
