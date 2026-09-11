#' Default calibration table shipped with rze
#'
#' A pre-computed AIC-threshold calibration table, used automatically by
#' \code{\link{rze}} when you do not supply your own. Generated with
#' \code{\link{calibrate_rze}} across a grid of tree size, sampling
#' fraction, and extinction fraction.
#'
#' @format A data frame with one row per calibration cell and columns:
#' \describe{
#'   \item{n}{Tree size (tip count).}
#'   \item{rho}{Sampling fraction.}
#'   \item{epsilon}{Extinction fraction.}
#'   \item{n_valid_replicates}{Null replicates that produced a usable value.}
#'   \item{aic_threshold}{The calibrated AIC-improvement acceptance threshold.}
#' }
#' @seealso \code{\link{calibrate_rze}}, \code{\link{resolve_calibration}}
"calibration_default"

#' rze: calibrated detection of clade-localized diversification shifts
#'
#' rze finds where on a time-calibrated phylogeny of extant taxa net
#' diversification rate shifts, and reports the net diversification rate of
#' each detected regime. It uses a stepwise search over candidate clades
#' with an AIC threshold calibrated against simulated constant-rate trees.
#'
#' @keywords internal
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics par plot plot.new plot.window rect axis title abline segments points text mtext layout image
#' @importFrom stats coef qnorm quantile reshape
#' @importFrom utils combn data
#' @importFrom tools file_ext
"_PACKAGE"
