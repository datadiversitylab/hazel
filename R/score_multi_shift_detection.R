# Score a detection result against multiple true shifts, extending the
# binary "found at least one" logic used in the single-shift validation
# into full per-shift accounting: how many of the N true shifts were
# recovered, and how many detected splits didn't correspond to any true
# shift (false positives within a tree that does have real shifts,
# distinct from the null-tree false positive rate, which is a different
# quantity, since here we're asking about precision, not just whether
# any signal exists).
#
# A true shift counts as recovered if some detected split's descendant
# clade overlaps that shift's true clade by at least overlap_threshold
# (proportion of the true clade's tips). A detected split counts as a
# true positive if it overlaps ANY true shift above that threshold, and
# as a false positive otherwise. This mirrors sensitivity/precision
# scoring generally, not something specific to this method.

score_multi_shift_detection <- function(result, tree, true_shifts,
                                         overlap_threshold = 0.7) {

  detected_clades <- lapply(result$splits, function(nd) {
    tree$tip.label[phangorn::Descendants(tree, nd, type = "tips")[[1]]]
  })

  # For each true shift, does any detected clade cover it well enough?
  true_shift_found <- vapply(true_shifts, function(shift) {
    true_tips <- shift$shifted_clade_tips
    if (length(detected_clades) == 0) return(FALSE)
    overlaps <- vapply(detected_clades, function(det_tips) {
      length(intersect(det_tips, true_tips)) / length(true_tips)
    }, numeric(1))
    any(overlaps >= overlap_threshold)
  }, logical(1))

  # For each detected split, does it correspond to any true shift, or is
  # it a false positive?
  detection_is_true_positive <- vapply(detected_clades, function(det_tips) {
    if (length(true_shifts) == 0) return(FALSE)
    overlaps <- vapply(true_shifts, function(shift) {
      true_tips <- shift$shifted_clade_tips
      length(intersect(det_tips, true_tips)) / length(true_tips)
    }, numeric(1))
    any(overlaps >= overlap_threshold)
  }, logical(1))

  n_true <- length(true_shifts)
  n_detected <- length(result$splits)
  n_recovered <- sum(true_shift_found)
  n_false_positive <- sum(!detection_is_true_positive)

  list(
    n_true = n_true,
    n_detected = n_detected,
    n_recovered = n_recovered,
    n_false_positive = n_false_positive,
    true_shift_found = true_shift_found,
    fraction_recovered = if (n_true > 0) n_recovered / n_true else NA_real_,
    all_recovered = n_true > 0 && n_recovered == n_true,
    exact_match = n_true > 0 && n_recovered == n_true && n_false_positive == 0
  )
}
