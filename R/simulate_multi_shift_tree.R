# Simulate a tree with MULTIPLE independent clade shifts grafted onto one
# background tree, generalizing simulate_clade_shift_tree() to more than
# one shift. Each shift is grafted at a distinct, randomly chosen
# background tip, using the same edge-splitting logic verified in
# simulate_clade_shift_tree() (receptor stem and donor height split the
# original tip edge in half, guaranteeing ultrametricity by construction,
# not by chance).
#
# shifts is a list of lists, each with n_shifted, lambda_ratio, and
# optionally epsilon (defaults to the background epsilon if omitted).

simulate_multi_shift_tree <- function(n_background, shifts,
                                       lambda_background = 1,
                                       epsilon = 0.2, rho = 1) {

  background <- simulate_null_tree(n_background, rho = rho, epsilon = epsilon)

  if (length(shifts) > n_background) {
    stop("Cannot graft more shifts than there are background tips available")
  }

  graft_tips <- sample(background$tip.label, length(shifts))

  tree <- background
  shift_records <- vector("list", length(shifts))

  for (i in seq_along(shifts)) {

    spec <- shifts[[i]]
    shift_epsilon <- if (is.null(spec$epsilon)) epsilon else spec$epsilon

    graft_tip_index <- which(tree$tip.label == graft_tips[i])
    tip_edge_index <- which(tree$edge[, 2] == graft_tip_index)
    original_tip_edge <- tree$edge.length[tip_edge_index]

    receptor_stem <- original_tip_edge / 2
    target_depth <- original_tip_edge / 2

    receptor <- tree
    receptor$tip.label[graft_tip_index] <- "NA"
    receptor$edge.length[tip_edge_index] <- receptor_stem

    shifted_lambda <- lambda_background * spec$lambda_ratio
    shifted_mu <- shift_epsilon * shifted_lambda

    shifted_clade <- TreeSim::sim.bd.taxa(
      n = spec$n_shifted, numbsim = 1, lambda = shifted_lambda,
      mu = shifted_mu, frac = rho, complete = FALSE, stochsampling = FALSE
    )[[1]]

    shifted_clade$tip.label <- paste0("shift", i, "_", shifted_clade$tip.label)

    shifted_height <- max(ape::node.depth.edgelength(shifted_clade))
    shifted_clade$edge.length <- shifted_clade$edge.length *
      (target_depth / shifted_height)
    shifted_clade$root.edge <- 0

    tree <- phytools::paste.tree(receptor, shifted_clade)

    shift_records[[i]] <- list(
      graft_tip = graft_tips[i],
      shifted_clade_tips = shifted_clade$tip.label,
      lambda_ratio = spec$lambda_ratio,
      n_shifted = spec$n_shifted
    )
  }

  if (!ape::is.ultrametric(tree)) {
    warning("Multi-shift grafted tree failed the ultrametric check, inspect before using")
  }

  list(tree = tree, shifts = shift_records)
}
