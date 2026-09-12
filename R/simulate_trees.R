# Simulate a constant-rate birth-death tree, the null model for
# calibration and validation. lambda fixed at 1, mu = epsilon * lambda,
# so a single epsilon sets the whole regime.

#' Simulate a constant-rate birth-death tree
#'
#' @param n Number of extant tips.
#' @param rho Sampling fraction.
#' @param epsilon Extinction fraction (mu = epsilon times lambda).
#' @return An ultrametric \code{phylo} object.
#' @export
simulate_null_tree <- function(n, rho = 1, epsilon = 0) {
  
  if (n < 4) stop("n must be at least 4")
  if (rho <= 0 || rho > 1) stop("rho must be in (0, 1]")
  if (epsilon < 0 || epsilon >= 1) stop("epsilon must be in [0, 1)")
  
  lambda <- 1
  mu <- epsilon * lambda
  
  trees <- TreeSim::sim.bd.taxa(
    n = n, numbsim = 1, lambda = lambda, mu = mu,
    frac = rho, complete = FALSE, stochsampling = FALSE
  )
  
  trees[[1]]
}

# Simulate a tree with a clade-localized rate shift, for validation.
#
# Built on phytools::paste.tree(), a function purpose-built for exactly
# this: grafting a donor clade onto a receptor tree at a designated
# attachment tip, guaranteeing ultrametricity because both trees are
# scaled to meet cleanly at the graft point. This replaces an earlier,
# hand-rolled approach using ape::bind.tree() and manual position
# arithmetic, which failed the ultrametric check twice and was abandoned
# rather than patched further, given how much trouble similar
# hand-verified-by-us-alone logic caused earlier in this project.
#
# paste.tree()'s documented requirements: the receptor tree needs the
# attachment tip labeled "NA", and the donor clade needs a root.edge
# (even if zero length).

#' Simulate a tree with one clade-localized rate shift
#'
#' @param n_background Background tree size.
#' @param n_shifted Tips in the shifted clade.
#' @param lambda_background Background speciation rate.
#' @param lambda_ratio Shifted clade rate as a multiple of background.
#' @param epsilon Extinction fraction.
#' @param rho Sampling fraction.
#' @return A list with the grafted tree and the true shifted-clade tips.
#' @export
simulate_clade_shift_tree <- function(n_background, n_shifted,
                                      lambda_background = 1,
                                      lambda_ratio = 3,
                                      epsilon = 0.2, rho = 1) {
  
  background <- simulate_null_tree(n_background, rho = rho, epsilon = epsilon)
  root_age <- max(ape::node.depth.edgelength(background))
  
  graft_tip_label <- sample(background$tip.label, 1)
  graft_tip_index <- which(background$tip.label == graft_tip_label)
  
  # paste.tree() replaces the "NA"-labeled tip with the donor clade's
  # root, it does not insert the donor partway along the existing edge.
  # That tip's own edge already carries it to the correct depth (the
  # present), so the donor clade must have height ZERO at the graft, or
  # its own internal structure will be added on top of an edge that
  # already reaches the present, overshooting total tree height. This
  # was confirmed by direct inspection after an earlier version of this
  # function (using the tip's edge length as the target height) produced
  # a non-ultrametric tree with shifted tips too deep by exactly the
  # donor clade's rescaled height.
  #
  # So instead: keep the receptor tip's edge exactly as simulated, and
  # rescale the shifted clade to have height equal to that SAME edge
  # length, then subtract that edge from the receptor tip before
  # grafting, so the two together sum to the original depth, not exceed
  # it.
  tip_edge_index <- which(background$edge[, 2] == graft_tip_index)
  original_tip_edge <- background$edge.length[tip_edge_index]
  
  # Split the original tip edge: half stays as the receptor's stem to the
  # graft point, half becomes the target height for the donor clade. This
  # guarantees receptor_stem + donor_height == original_tip_edge exactly.
  receptor_stem <- original_tip_edge / 2
  target_depth <- original_tip_edge / 2
  
  receptor <- background
  receptor$tip.label[graft_tip_index] <- "NA"
  receptor$edge.length[tip_edge_index] <- receptor_stem
  
  shifted_lambda <- lambda_background * lambda_ratio
  shifted_mu <- epsilon * shifted_lambda
  
  shifted_clade <- TreeSim::sim.bd.taxa(
    n = n_shifted, numbsim = 1, lambda = shifted_lambda, mu = shifted_mu,
    frac = rho, complete = FALSE, stochsampling = FALSE
  )[[1]]
  
  # TreeSim assigns generic tip labels (t1, t2, ...) independently to
  # each simulated tree, so background and shifted_clade will very likely
  # collide on names. Relabel the donor clade's tips to guarantee
  # uniqueness before grafting.
  shifted_clade$tip.label <- paste0("shift_", shifted_clade$tip.label)
  
  shifted_height <- max(ape::node.depth.edgelength(shifted_clade))
  shifted_clade$edge.length <- shifted_clade$edge.length *
    (target_depth / shifted_height)
  shifted_clade$root.edge <- 0
  
  grafted <- phytools::paste.tree(receptor, shifted_clade)
  
  if (!ape::is.ultrametric(grafted)) {
    warning("Grafted tree failed the ultrametric check, inspect before using")
  }
  
  list(tree = grafted, graft_tip = graft_tip_label,
       shifted_clade_tips = shifted_clade$tip.label,
       attach_depth = root_age - receptor_stem)
}
#' Simulate a tree with a nested rate shift
#'
#' Builds a tree with one rate shift nested inside another: an outer clade
#' with its own diversification rate that itself contains an inner clade
#' with a further rate change. This is the harder case for shift detection,
#' since the inner and outer regimes must be attributed correctly at once.
#'
#' The inner shift is grafted onto a tip of the outer shifted clade, so the
#' nesting is real by construction. The two clades' tips are labeled so the
#' inner tips are a strict subset of the outer clade's descendants.
#'
#' @param n_background Background tree size (tip count).
#' @param n_outer Tips in the outer shifted clade before the inner graft.
#' @param n_inner Tips in the inner (nested) shifted clade.
#' @param lambda_background Background speciation rate.
#' @param outer_ratio Outer clade speciation rate as a multiple of background.
#' @param inner_ratio Inner clade speciation rate as a multiple of background.
#' @param epsilon Extinction fraction.
#' @param rho Sampling fraction.
#' @return A list with the tree, the outer clade's tip set, and the inner
#'   clade's tip set (a subset of the outer set).
#' @export
simulate_nested_shift_tree <- function(n_background, n_outer, n_inner,
                                       lambda_background = 1,
                                       outer_ratio = 3, inner_ratio = 6,
                                       epsilon = 0.2, rho = 1) {

  # Graft the outer shift onto a background tip, using the same
  # edge-splitting graft as the single-shift simulator.
  graft_shift <- function(receptor_tree, target_tip, n_shifted,
                          lambda_ratio, prefix) {
    idx <- which(receptor_tree$tip.label == target_tip)
    edge_i <- which(receptor_tree$edge[, 2] == idx)
    orig_edge <- receptor_tree$edge.length[edge_i]

    receptor_stem <- orig_edge / 2
    target_depth <- orig_edge / 2

    receptor <- receptor_tree
    receptor$tip.label[idx] <- "NA"
    receptor$edge.length[edge_i] <- receptor_stem

    lam <- lambda_background * lambda_ratio
    mu <- epsilon * lam
    clade <- TreeSim::sim.bd.taxa(
      n = n_shifted, numbsim = 1, lambda = lam, mu = mu,
      frac = rho, complete = FALSE, stochsampling = FALSE
    )[[1]]
    clade$tip.label <- paste0(prefix, clade$tip.label)

    h <- max(ape::node.depth.edgelength(clade))
    clade$edge.length <- clade$edge.length * (target_depth / h)
    clade$root.edge <- 0

    list(tree = phytools::paste.tree(receptor, clade),
         clade_tips = clade$tip.label)
  }

  background <- simulate_null_tree(n_background, rho = rho, epsilon = epsilon)

  # Outer shift onto a random background tip
  outer_target <- sample(background$tip.label, 1)
  outer <- graft_shift(background, outer_target, n_outer, outer_ratio, "outer_")
  tree <- outer$tree
  outer_tips <- outer$clade_tips

  # Inner shift onto a tip WITHIN the outer clade, so it nests
  inner_target <- sample(outer_tips, 1)
  inner <- graft_shift(tree, inner_target, n_inner, inner_ratio, "inner_")
  tree <- inner$tree
  inner_tips <- inner$clade_tips

  # The outer clade's descendant tips now include the inner tips plus the
  # remaining outer tips (all but the one replaced by the inner graft).
  outer_tips_final <- c(setdiff(outer_tips, inner_target), inner_tips)

  if (!ape::is.ultrametric(tree)) {
    warning("Nested grafted tree failed the ultrametric check, inspect before using")
  }

  list(tree = tree,
       outer_clade_tips = outer_tips_final,
       inner_clade_tips = inner_tips,
       outer_ratio = outer_ratio,
       inner_ratio = inner_ratio)
}
