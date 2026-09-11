# Extract per-regime net diversification rates from a completed clade
# search.
#
# The package reports net diversification (lambda minus mu) as the primary
# rate, because it is far better identified from an extant-only tree than
# lambda and mu separately. The underlying fit still estimates lambda and
# mu (that is how make.bd.split works), but on small or data-poor clades
# the optimizer routinely drives mu to its lower boundary, at which point
# the separation of speciation from extinction is not supported by the
# data even though net diversification usually still is. The boundary_mu
# flag travels with every regime so the user always knows when this
# happened, regardless of how they choose to see the components.
#
# on_boundary controls what the user sees for a regime whose mu hit the
# boundary:
#   "flag"   report the value, mark it (default)
#   "na"     report NA for that regime's rate, since the fit is unreliable
#   "keep"   report the value with no special treatment

regime_rates <- function(result, tree, on_boundary = c("flag", "na", "keep")) {

  on_boundary <- match.arg(on_boundary)

  final <- result$fit
  splits <- result$splits

  # Map each partition to the clade it represents. make.bd.split orders
  # partitions as: partition 1 is the root (background) regime, partition
  # i+1 corresponds to nodes[i], where nodes is what was passed to
  # make.bd.split. IMPORTANT: fit_split_tree passes result$splits already
  # sorted, and the final fit stores that same sorted node vector in
  # final$nodes, so the ordering here is anchored to final$nodes rather
  # than to result$splits directly. This mapping is pinned by a package
  # test (test-partition-order); if diversitree ever changes its ordering
  # convention, that test fails loudly rather than silently mislabeling
  # which clade got which rate.
  ordered_nodes <- if (!is.null(final$nodes)) final$nodes else splits
  splits <- ordered_nodes
  n_regimes <- length(final$lambda)

  net_div <- final$lambda - final$mu
  boundary <- final$mu < 1e-6

  regime_label <- c("background", if (length(splits) > 0) paste0("shift_", seq_along(splits)))
  regime_node <- c(NA_integer_, splits)

  clade_size <- vapply(seq_len(n_regimes), function(i) {
    if (i == 1) {
      # background regime size is total tips minus all shifted clades, an
      # approximation, since split regimes carve tips away from the root
      length(tree$tip.label) - sum(vapply(splits, function(nd) {
        length(phangorn::Descendants(tree, nd, type = "tips")[[1]])
      }, integer(1)))
    } else {
      length(phangorn::Descendants(tree, splits[i - 1], type = "tips")[[1]])
    }
  }, integer(1))

  reported_rate <- net_div
  if (on_boundary == "na") {
    reported_rate[boundary] <- NA_real_
  }

  out <- data.frame(
    regime = regime_label,
    node = regime_node,
    n_tips = clade_size,
    net_diversification = reported_rate,
    boundary_mu = boundary,
    stringsAsFactors = FALSE
  )

  attr(out, "on_boundary") <- on_boundary
  attr(out, "lambda") <- final$lambda
  attr(out, "mu") <- final$mu
  out
}
