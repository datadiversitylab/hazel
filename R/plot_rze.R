# The main rze plot: a tree with every branch colored by the net
# diversification rate of the regime it belongs to, detected shifts
# marked at their nodes, and a rate colorbar beneath. This is the figure
# meant to communicate the whole result at a glance.
#
# Branch-to-regime assignment: a branch belongs to the most recent
# (deepest-nested) shift regime that contains it. Concretely, for each
# edge, we find which detected shift clades contain that edge's child
# node, and assign the smallest such clade; edges in no shift clade get
# the background regime.

plot.rze_result <- function(x, show_tip_labels = FALSE,
                            shift_cex = 2, edge_width = 3,
                            layout = TRUE, ...) {

  tree <- x$tree
  shifts <- x$shifts
  rates <- x$rates

  # Rate per regime, indexed the same way rates rows are: row 1 is
  # background, rows 2+ are shift regimes in order of shifts
  regime_rate <- rates$net_diversification
  names(regime_rate) <- rates$regime

  # For each edge, determine its regime. Start everything as background,
  # then for each shift clade (processed largest first so nested smaller
  # clades overwrite), color edges within that clade.
  n_edges <- nrow(tree$edge)
  edge_regime <- rep(1L, n_edges)  # 1 = background row in rates

  if (length(shifts) > 0) {
    clade_tip_sets <- lapply(shifts, function(nd) {
      phangorn::Descendants(tree, nd, type = "tips")[[1]]
    })
    clade_order <- order(vapply(clade_tip_sets, length, integer(1)),
                          decreasing = TRUE)

    # A node belongs to shift i if it is that shift's node or a descendant
    descendant_nodes <- function(tree, node) {
      out <- integer(0)
      to_visit <- node
      while (length(to_visit) > 0) {
        cur <- to_visit[1]; to_visit <- to_visit[-1]
        out <- c(out, cur)
        kids <- tree$edge[tree$edge[, 1] == cur, 2]
        to_visit <- c(to_visit, kids)
      }
      out
    }

    for (i in clade_order) {
      nodes_in_clade <- descendant_nodes(tree, shifts[i])
      in_clade_edge <- tree$edge[, 2] %in% nodes_in_clade
      edge_regime[in_clade_edge] <- i + 1L  # +1 because row 1 is background
    }
  }

  # Color scale spans all regime rates
  zlim <- range(regime_rate, na.rm = TRUE)
  pal <- rze_palette(256)
  edge_colors <- rze_map_color(regime_rate[edge_regime], palette = pal, zlim = zlim)

  # Layout: tree on top, colorbar beneath
  if (layout) {
    old_layout <- graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(6, 1))
    on.exit(graphics::layout(1), add = TRUE)
  }

  old_par <- graphics::par(mar = c(0.5, 0.5, 2.5, 0.5), bg = "white")
  on.exit(graphics::par(old_par), add = TRUE)

  ape::plot.phylo(
    tree,
    edge.color = edge_colors,
    edge.width = edge_width,
    show.tip.label = show_tip_labels,
    cex = 0.5,
    no.margin = FALSE,
    main = ""
  )
  graphics::title(main = sprintf("rze: %d diversification shift%s detected",
                                 length(shifts),
                                 if (length(shifts) == 1) "" else "s"),
                  col.main = "#111111", font.main = 1, cex.main = 1.3)

  # Mark shifts
  if (length(shifts) > 0) {
    ape::nodelabels(node = shifts, pch = 21, cex = shift_cex,
                     bg = rze_shift_color(), col = "#003030", lwd = 2)
  }

  if (layout) {
    rze_colorbar(zlim, palette = pal, label = "net diversification rate")
  }

  invisible(x)
}
