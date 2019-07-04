#' Simulated data set for use in examples
#'
#' A dataset containing a simulated response variable and several predictor
#'   variables for use in demonstrating occupancy-detection models.
#'
#' @format A data frame with 200 rows and 10 variables:
#' \describe{
#'   \item{response}{binary detection data at a site and survey}
#'   \item{occ_predictor1}{simulated predictor variable that affects occupancy, recorded at site level}
#'   \item{occ_predictor2}{simulated predictor variable that affects occupancy, recorded at site level}
#'   \item{occ_random1}{simulated grouping variable that affects occupancy, recorded at site level}
#'   \item{occ_random2}{simulated grouping variable that affects occupancy, recorded at site level}
#'   \item{detect_predictor1}{simulated predictor variable that affects detection, recorded at survey level}
#'   \item{detect_predictor2}{simulated predictor variable that affects detection, recorded at survey level}
#'   \item{detect_random1}{simulated grouping variable that affects detection, recorded at survey level}
#'   \item{site}{site identifier}
#'   \item{survey}{survey identifier (visit number to a given site)}
#' }
"occupancy_data"