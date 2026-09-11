# Forward-then-backward stepwise search for clade-localized rate shifts,
# using fit_split_tree() (diversitree::make.bd.split) as the engine.
#
# Two modes, controlled by aic_threshold, same convention as before.
#
# Normal mode: aic_threshold is a positive number, the calibrated minimum
# AIC improvement required to accept a split node, from
# calibrate_medusa_null(). Runs forward until no candidate clears the
# threshold, then prunes backward.
#
# Calibration mode: aic_threshold is NULL. Runs a single forward step
# only, returns the best candidate and its AIC improvement without
# accepting it. This is what calibrate_medusa_null() calls repeatedly to
# build the null distribution the threshold is calibrated against.

#' Stepwise forward-backward clade search
#'
#' Adds shift regimes one at a time while AIC improves past the threshold,
#' then prunes. With \code{aic_threshold = NULL} it runs a single
#' calibration step instead, used by \code{\link{calibrate_rze}}.
#'
#' @param tree A \code{phylo} object.
#' @param candidates Candidate node indices.
#' @param aic_threshold Acceptance threshold, or \code{NULL} for
#'   calibration mode.
#' @param max_splits Cap on accepted shifts.
#' @param sampling.f Sampling fraction.
#' @param verbose Whether to print progress.
#' @param ... Passed to the underlying fit.
#' @return In normal mode, a list with accepted splits, the fit, and
#'   history; in calibration mode, the best candidate and its improvement.
#' @export
stepwise_clade_search <- function(tree, candidates, aic_threshold = NULL,
                                   max_splits = NULL, sampling.f = 1,
                                   verbose = TRUE) {

  if (verbose) cat("Fitting whole-tree baseline...\n")
  baseline <- fit_whole_tree(tree, sampling.f = sampling.f)
  baseline_aic <- baseline$AIC
  if (verbose) cat(sprintf("  baseline AIC = %.2f (lambda = %.3f, mu = %.3f)\n",
                            baseline_aic, baseline$lambda, baseline$mu))

  fit_at <- function(nodes) {
    tryCatch(
      fit_split_tree(tree, nodes, sampling.f = sampling.f),
      error = function(e) list(AIC = Inf, boundary_mu = NA, error = conditionMessage(e))
    )
  }

  # Calibration mode: one forward step, report best candidate, accept
  # nothing
  if (is.null(aic_threshold)) {

    if (length(candidates) == 0) {
      return(list(best_node = NA_integer_, aic_improvement = 0))
    }

    if (verbose) cat(sprintf("Scanning %d candidates (calibration step)...\n",
                              length(candidates)))

    results <- lapply(seq_along(candidates), function(i) {
      if (verbose) cat(sprintf("  [%d/%d] node = %d\n", i, length(candidates),
                                candidates[i]))
      fit_at(candidates[i])
    })

    aics <- vapply(results, function(r) r$AIC, numeric(1))
    best_index <- which.min(aics)

    return(list(
      best_node = candidates[best_index],
      aic_improvement = baseline_aic - aics[best_index],
      boundary_mu = results[[best_index]]$boundary_mu
    ))
  }

  # Normal mode: forward search
  accepted <- integer(0)
  remaining <- candidates
  history <- data.frame(node = integer(0), aic = numeric(0),
                         improvement = numeric(0), boundary_mu = logical(0))

  if (is.null(max_splits)) {
    max_splits <- max(1, floor(length(tree$tip.label) / 20))
  }

  step_number <- 0
  while (length(remaining) > 0 && length(accepted) < max_splits) {

    step_number <- step_number + 1
    if (verbose) {
      cat(sprintf("Forward step %d: scanning %d candidates, %d already accepted...\n",
                   step_number, length(remaining), length(accepted)))
    }

    results <- lapply(seq_along(remaining), function(i) {
      if (verbose) cat(sprintf("  [%d/%d] node = %d\n", i, length(remaining),
                                remaining[i]))
      fit_at(c(accepted, remaining[i]))
    })

    aics <- vapply(results, function(r) r$AIC, numeric(1))
    best_index <- which.min(aics)
    improvement <- baseline_aic - aics[best_index]

    if (verbose) {
      cat(sprintf("  best candidate: node = %d, AIC improvement = %.2f\n",
                   remaining[best_index], improvement))
    }

    if (improvement < aic_threshold) {
      if (verbose) cat("  below threshold, stopping forward search\n")
      break
    }

    accepted <- c(accepted, remaining[best_index])
    baseline_aic <- aics[best_index]
    history <- rbind(history, data.frame(
      node = remaining[best_index], aic = baseline_aic,
      improvement = improvement,
      boundary_mu = isTRUE(results[[best_index]]$boundary_mu)
    ))
    remaining <- remaining[-best_index]
    if (verbose) cat(sprintf("  accepted split at node %d\n", accepted[length(accepted)]))
  }

  # Backward pruning
  if (verbose) cat("Backward pruning...\n")
  repeat {
    if (length(accepted) == 0) break

    drop_results <- lapply(accepted, function(nd) fit_at(setdiff(accepted, nd)))
    drop_aic <- vapply(drop_results, function(r) r$AIC, numeric(1))
    worst_loss <- drop_aic - baseline_aic
    drop_index <- which.min(worst_loss)

    if (worst_loss[drop_index] > aic_threshold) break

    if (verbose) cat(sprintf("  dropping split at node %d\n", accepted[drop_index]))
    accepted <- accepted[-drop_index]
    baseline_aic <- drop_aic[drop_index]
  }
  if (verbose) cat("Done.\n")

  final_fit <- if (length(accepted) > 0) fit_at(accepted) else baseline

  list(
    splits = sort(accepted),
    fit = final_fit,
    history = history
  )
}
