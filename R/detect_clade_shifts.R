# User-facing entry point. Runs a calibrated stepwise clade search on a
# real tree, using a calibration table produced by calibrate_medusa_null()
# (or the shipped default, once one exists).

detect_clade_shifts <- function(tree, calibration_table, rho = 1,
                                 epsilon = 0.2, min_clade_size = 5,
                                 verbose = TRUE) {

  n <- length(tree$tip.label)

  threshold <- get_calibrated_threshold(calibration_table, n, rho, epsilon)
  if (verbose) cat(sprintf("Using calibrated AIC threshold: %.2f\n", threshold))

  candidates <- candidate_nodes(tree, min_clade_size = min_clade_size)
  if (verbose) cat(sprintf("%d candidate split nodes (min clade size %d)\n",
                            length(candidates), min_clade_size))

  result <- stepwise_clade_search(
    tree, candidates, aic_threshold = threshold,
    sampling.f = rho, verbose = verbose
  )

  result$calibration_used <- list(n = n, rho = rho, epsilon = epsilon,
                                   threshold = threshold)
  result
}
