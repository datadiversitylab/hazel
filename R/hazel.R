#' Detect diversification rate shifts and estimate per-regime rates
#'
#' Runs the full hazel pipeline on a tree: calibrated clade-shift detection,
#' then per-regime net diversification rate estimation. Returns a classed
#' \code{hazel_result} object with print, summary, and plot methods.
#'
#' The calibration table is resolved through \code{\link{resolve_calibration}}:
#' a table you pass in wins, otherwise the table shipped with the package is
#' used, and only if neither is available does hazel calibrate on the fly for
#' your specific tree (which is slower).
#'
#' hazel reports net diversification (speciation minus extinction) per regime.
#' On an extant-only tree, speciation and extinction are not separately
#' identifiable in general, so net diversification is the honest primary
#' output. When a regime's extinction estimate hits its lower boundary, that
#' regime is flagged, since the speciation/extinction split is not supported
#' by the data there even though net diversification usually still is.
#'
#' @param tree An ultrametric \code{phylo} object of extant taxa.
#' @param calibration A calibration table (a data frame from
#'   \code{\link{calibrate_hazel}}), a path to a saved \code{.rds}/\code{.rda}
#'   file, or \code{NULL} to use the shipped default (or on-the-fly
#'   calibration if no default is available).
#' @param rho Sampling fraction, the proportion of species included in the
#'   tree. A single global value.
#' @param epsilon Assumed extinction fraction, used to pick the calibration
#'   cell.
#' @param min_clade_size Smallest clade the search may propose as its own
#'   regime.
#' @param on_boundary How to report a regime whose extinction estimate hit
#'   the boundary: \code{"flag"} reports the value and marks it,
#'   \code{"na"} returns \code{NA} for that regime, \code{"keep"} reports it
#'   with no special treatment.
#' @param verbose Whether to print progress messages.
#'
#' @return An \code{hazel_result} object: a list with the tree, detected shift
#'   nodes, a per-regime rate table, the search history, the calibration
#'   used, the underlying fit, and the settings.
#'
#' @examples
#' \dontrun{
#' # With the shipped calibration table (no calibration step needed)
#' result <- hazel(tree, rho = 0.8)
#' print(result)
#' plot(result)
#' plot_regime_rates(result)
#' }
#'
#' @seealso \code{\link{calibrate_hazel}}, \code{\link{resolve_calibration}},
#'   \code{\link{plot_regime_rates}}
#' @export
hazel <- function(tree, calibration = NULL, rho = 1, epsilon = 0.2,
                min_clade_size = 5, on_boundary = c("flag", "na", "keep"),
                verbose = TRUE) {

  on_boundary <- match.arg(on_boundary)

  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  if (!ape::is.ultrametric(tree)) {
    stop("tree must be ultrametric (a time-calibrated tree of extant taxa)")
  }

  if (is.null(calibration) || is.character(calibration) || is.data.frame(calibration)) {
    resolved <- resolve_calibration(calibration, verbose = verbose)
    if (is.null(resolved)) {
      if (verbose) message("No calibration supplied or shipped, calibrating on the fly (this is slower)...")
      resolved <- calibrate_hazel(
        n_grid = length(tree$tip.label),
        rho_grid = rho, epsilon_grid = epsilon,
        n_replicates = 50, min_clade_size = min_clade_size, verbose = FALSE
      )
    }
    calibration <- resolved
  }

  detection <- detect_clade_shifts(
    tree, calibration, rho = rho, epsilon = epsilon,
    min_clade_size = min_clade_size, verbose = verbose
  )

  rates <- regime_rates(detection, tree, on_boundary = on_boundary)

  structure(
    list(
      tree = tree,
      shifts = detection$splits,
      rates = rates,
      history = detection$history,
      calibration_used = detection$calibration_used,
      fit = detection$fit,
      settings = list(rho = rho, epsilon = epsilon,
                       min_clade_size = min_clade_size,
                       on_boundary = on_boundary)
    ),
    class = "hazel_result"
  )
}

#' @param x An \code{hazel_result} object.
#' @param ... Further arguments (unused).
#' @rdname hazel
#' @export
print.hazel_result <- function(x, ...) {
  n_shifts <- length(x$shifts)
  cat("hazel diversification shift analysis\n")
  cat(sprintf("  tree: %d tips\n", length(x$tree$tip.label)))
  cat(sprintf("  sampling fraction (rho): %.2f\n", x$settings$rho))
  cat(sprintf("  detected shifts: %d\n", n_shifts))
  if (n_shifts > 0) {
    cat("\n")
    rates_show <- x$rates
    rates_show$net_diversification <- round(rates_show$net_diversification, 4)
    print(rates_show, row.names = FALSE)
    if (any(x$rates$boundary_mu)) {
      cat("\n  note: regimes flagged boundary_mu had extinction pinned at the\n")
      cat("  boundary; net diversification is still reported but the\n")
      cat("  speciation/extinction split is not supported by the data there.\n")
    }
  } else {
    cat("  no rate heterogeneity detected; whole-tree rate applies\n")
  }
  invisible(x)
}

#' @param object An \code{hazel_result} object.
#' @rdname hazel
#' @export
summary.hazel_result <- function(object, ...) {
  print(object, ...)
}
