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

# Interpolate a calibrated threshold for a tree's actual n, rho, epsilon
# against a calibration table, nearest-neighbor for now, simplest thing
# that works, refine to real interpolation once the shipped grid exists.

get_calibrated_threshold <- function(calibration_table, n, rho, epsilon) {

  d <- sqrt(
    ((calibration_table$n - n) / max(calibration_table$n))^2 +
    ((calibration_table$rho - rho) / max(calibration_table$rho, 1))^2 +
    ((calibration_table$epsilon - epsilon) / max(calibration_table$epsilon, 1))^2
  )

  calibration_table$aic_threshold[which.min(d)]
}

# Public-facing name for the calibration function. calibrate_medusa_null()
# remains as the internal name the rest of the code and tests already use;
# calibrate_rze() is the name users see and call.
calibrate_rze <- function(...) calibrate_medusa_null(...)
