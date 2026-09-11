# Main user-facing entry point for rze.
#
# Runs the full pipeline on a tree: calibrated clade-shift detection, then
# per-regime net diversification rate estimation. Returns a classed
# rze_result object with print, summary, and plot methods.
#
# tree            an ultrametric phylo object of extant taxa
# calibration     a calibration table from calibrate_rze() (or the shipped
#                 default); if NULL, calibration is run on the fly (slower)
# rho             sampling fraction, a single value (global) for now
# epsilon         assumed extinction fraction for calibration cell lookup
# min_clade_size  smallest clade the search may propose as its own regime
# on_boundary     how to report a regime whose extinction estimate hit the
#                 boundary, passed to regime_rates()

rze <- function(tree, calibration = NULL, rho = 1, epsilon = 0.2,
                min_clade_size = 5, on_boundary = c("flag", "na", "keep"),
                verbose = TRUE) {

  on_boundary <- match.arg(on_boundary)

  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  if (!ape::is.ultrametric(tree)) {
    stop("tree must be ultrametric (a time-calibrated tree of extant taxa)")
  }

  if (is.null(calibration)) {
    if (verbose) message("No calibration supplied, calibrating on the fly (this is slower)...")
    calibration <- calibrate_rze(
      n_grid = length(tree$tip.label),
      rho_grid = rho, epsilon_grid = epsilon,
      n_replicates = 50, min_clade_size = min_clade_size, verbose = FALSE
    )
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
    class = "rze_result"
  )
}

print.rze_result <- function(x, ...) {
  n_shifts <- length(x$shifts)
  cat("rze diversification shift analysis\n")
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

summary.rze_result <- function(object, ...) {
  print(object, ...)
}
