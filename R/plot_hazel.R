#' Plot a hazel result as a branch-colored tree
#'
#' A phylogeny with branches colored by regime.
#'
#' @param x A \code{hazel_result} object.
#' @param show_tip_labels Whether to draw tip labels.
#' @param mark_size Size of the shift markers.
#' @param edge_width Branch line width.
#' @param key Whether to draw the per-regime rate key.
#' @param ... Unused.
#' @return Called for its plot side effect; returns \code{x} invisibly.
#' @export
plot.hazel_result <- function(x, show_tip_labels = FALSE,
                              mark_size = 1.8, edge_width = 2.6,
                              key = TRUE, ...) {
  
  tree <- x$tree
  shifts <- x$shifts
  rates <- x$rates
  
  n_regimes <- nrow(rates)
  
  # Distinct color per regime, recycled from the hazel regime palette if
  # there are more regimes than colors. Regime 1 is the background.
  regime_pal <- hazel_regime_colors()
  regime_col <- regime_pal[((seq_len(n_regimes) - 1) %% length(regime_pal)) + 1]
  
  # Assign each edge to its regime. Everything starts as background; each
  # shift clade (largest first, so nested clades overwrite) claims its edges.
  n_edges <- nrow(tree$edge)
  edge_regime <- rep(1L, n_edges)
  
  if (length(shifts) > 0) {
    clade_tip_sets <- lapply(shifts, function(nd) {
      phangorn::Descendants(tree, nd, type = "tips")[[1]]
    })
    clade_order <- order(vapply(clade_tip_sets, length, integer(1)),
                         decreasing = TRUE)
    
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
      edge_regime[in_clade_edge] <- i + 1L
    }
  }
  
  edge_colors <- regime_col[edge_regime]
  
  old_par <- graphics::par(mar = c(1, 1, 1, 1), bg = "white")
  on.exit(graphics::par(old_par), add = TRUE)
  
  ape::plot.phylo(
    tree,
    edge.color = edge_colors,
    edge.width = edge_width,
    show.tip.label = show_tip_labels,
    cex = 0.5,
    main = ""
  )
  
  if (length(shifts) > 0) {
    ape::nodelabels(node = shifts, pch = 21, cex = mark_size,
                    bg = hazel_mark(), col = hazel_mark_border(), lwd = 1.8)
  }
  
  if (key) {
    hazel_regime_key(rates$regime, rates$net_diversification, regime_col,
                     rates$boundary_mu)
  }
  
  invisible(x)
}
