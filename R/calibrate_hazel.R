# Calibrate the AIC-improvement threshold against a target false-positive
# rate, by simulating null (no true shift) trees and running one
# uncalibrated forward step of stepwise_clade_search() on each.
#
# Returns one row per (n, rho, epsilon) cell: the target-percentile AIC
# improvement under the null, to be used as aic_threshold in normal-mode
# stepwise_clade_search() calls on real trees matching that cell.

calibrate_medusa_null <- function(n_grid, rho_grid, epsilon_grid,
                                  n_replicates, min_clade_size = 5,
                                  target_percentile = 0.95, verbose = TRUE) {
  
  grid <- expand.grid(n = n_grid, rho = rho_grid, epsilon = epsilon_grid)
  results <- vector("list", nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    
    n <- grid$n[i]; rho <- grid$rho[i]; epsilon <- grid$epsilon[i]
    if (verbose) {
      cat(sprintf("\nCell %d/%d: n=%d, rho=%.2f, epsilon=%.2f\n",
                  i, nrow(grid), n, rho, epsilon))
    }
    
    improvements <- numeric(n_replicates)
    
    for (r in seq_len(n_replicates)) {
      tree <- simulate_null_tree(n = n, rho = rho, epsilon = epsilon)
      candidates <- candidate_nodes(tree, min_clade_size = min_clade_size)
      
      if (length(candidates) == 0) {
        improvements[r] <- NA
        next
      }
      
      step <- tryCatch(
        stepwise_clade_search(tree, candidates, aic_threshold = NULL,
                              sampling.f = rho, verbose = FALSE),
        error = function(e) list(aic_improvement = NA)
      )
      improvements[r] <- step$aic_improvement
      
      if (verbose && r %% max(1, floor(n_replicates / 10)) == 0) {
        cat(sprintf("  replicate %d/%d\n", r, n_replicates))
      }
    }
    
    valid <- improvements[!is.na(improvements)]
    threshold <- if (length(valid) > 0) {
      stats::quantile(valid, probs = target_percentile, na.rm = TRUE)
    } else {
      NA
    }
    
    results[[i]] <- data.frame(
      n = n, rho = rho, epsilon = epsilon,
      n_valid_replicates = length(valid),
      aic_threshold = unname(threshold)
    )
  }
  
  do.call(rbind, results)
}

#' Look up the calibrated AIC threshold for a tree
#'
#' Given a calibration table and a tree's size, sampling fraction, and
#' extinction fraction, returns the AIC-improvement threshold for the
#' nearest matching cell. Used internally by \code{\link{hazel}}, and exported
#' for custom calibration workflows.
#'
#' @param calibration_table A calibration table from \code{\link{calibrate_hazel}}.
#' @param n Tree size (tip count).
#' @param rho Sampling fraction.
#' @param epsilon Extinction fraction.
#' @return The calibrated AIC-improvement threshold for the nearest cell.
#' @export
get_calibrated_threshold <- function(calibration_table, n, rho, epsilon) {
  
  d <- sqrt(
    ((calibration_table$n - n) / max(calibration_table$n))^2 +
      ((calibration_table$rho - rho) / max(calibration_table$rho, 1))^2 +
      ((calibration_table$epsilon - epsilon) / max(calibration_table$epsilon, 1))^2
  )
  
  calibration_table$aic_threshold[which.min(d)]
}

#' Calibrate the AIC threshold against simulated null trees
#'
#' Simulates constant-rate (no-shift) trees across a grid of tree size,
#' sampling fraction, and extinction fraction, and records the AIC
#' improvement a spurious shift achieves under the null. The chosen
#' percentile of that null distribution becomes the acceptance threshold
#' for each cell, so the false-positive rate is controlled rather than left
#' to an uncalibrated default.
#'
#' @param n_grid Tree sizes (tip counts) to calibrate for.
#' @param rho_grid Sampling fractions to calibrate for.
#' @param epsilon_grid Extinction fractions to calibrate for.
#' @param n_replicates Null trees simulated per cell. More is more reliable;
#'   small counts give noisy thresholds.
#' @param min_clade_size Smallest clade the search may propose.
#' @param target_percentile Null percentile used as the threshold
#'   (0.95 targets a 5 percent false-positive rate).
#' @param verbose Whether to print progress.
#'
#' @return A data frame with one row per cell: \code{n}, \code{rho},
#'   \code{epsilon}, the number of valid replicates, and
#'   \code{aic_threshold}.
#' @export
calibrate_hazel <- function(n_grid, rho_grid, epsilon_grid,
                          n_replicates, min_clade_size = 5,
                          target_percentile = 0.95, verbose = TRUE) {
  calibrate_medusa_null(n_grid = n_grid, rho_grid = rho_grid,
                        epsilon_grid = epsilon_grid,
                        n_replicates = n_replicates,
                        min_clade_size = min_clade_size,
                        target_percentile = target_percentile,
                        verbose = verbose)
}
