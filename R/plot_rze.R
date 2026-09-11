#' Plot an rze result as a branch-colored tree
#'
#' The main rze figure: a phylogeny with every branch colored by the net
#' diversification rate of the regime it belongs to, and detected shifts
#' marked at their nodes. A compact rate legend is drawn inside the plot.
#' Built in base R, single panel, so it is robust across devices.
#'
#' @param x An \code{rze_result} object.
#' @param show_tip_labels Whether to draw tip labels.
#' @param shift_cex Size of the shift markers.
#' @param edge_width Branch line width.
#' @param legend Whether to draw the inline rate legend.
#' @param ... Unused.
#' @return Called for its plot side effect; returns \code{x} invisibly.
#' @export
plot.rze_result <- function(x, show_tip_labels = FALSE,
                            shift_cex = 2, edge_width = 3,
                            legend = TRUE, ...) {

  tree <- x$tree
  shifts <- x$shifts
  rates <- x$rates

  # Rate per regime: row 1 is background, rows 2+ are shift regimes in the
  # order of shifts
  regime_rate <- rates$net_diversification
  names(regime_rate) <- rates$regime

  # Assign each edge to its regime. Everything starts as background; each
  # shift clade (largest first, so nested clades overwrite) colors its edges.
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

  zlim <- range(regime_rate, na.rm = TRUE)
  pal <- rze_palette(256)
  edge_colors <- rze_map_color(regime_rate[edge_regime], palette = pal, zlim = zlim)

  # Single panel. Save and restore only margins, no layout, so there is
  # nothing to corrupt the graphics state.
  old_par <- graphics::par(mar = c(1, 1, 3, 1), bg = "white")
  on.exit(graphics::par(old_par), add = TRUE)

  ape::plot.phylo(
    tree,
    edge.color = edge_colors,
    edge.width = edge_width,
    show.tip.label = show_tip_labels,
    cex = 0.5,
    main = ""
  )
  graphics::title(main = sprintf("rze: %d diversification shift%s detected",
                                 length(shifts),
                                 if (length(shifts) == 1) "" else "s"),
                  col.main = "#111111", font.main = 1, cex.main = 1.3)

  if (length(shifts) > 0) {
    ape::nodelabels(node = shifts, pch = 21, cex = shift_cex,
                    bg = rze_shift_color(), col = "#003030", lwd = 2)
  }

  # Inline legend: a small color ramp with low/high labels, drawn in the
  # top-left in user coordinates. No extra panel, no layout.
  if (legend && is.finite(zlim[1]) && is.finite(zlim[2])) {
    usr <- graphics::par("usr")
    xr <- usr[2] - usr[1]
    yr <- usr[4] - usr[3]

    bx0 <- usr[1] + 0.02 * xr
    bx1 <- usr[1] + 0.22 * xr
    by1 <- usr[4] - 0.02 * yr
    by0 <- usr[4] - 0.05 * yr

    n <- length(pal)
    xs <- seq(bx0, bx1, length.out = n + 1)
    graphics::rect(xs[-(n + 1)], by0, xs[-1], by1, col = pal, border = NA)
    graphics::rect(bx0, by0, bx1, by1, border = rze_neutral(), lwd = 1)

    graphics::text(bx0, by0 - 0.015 * yr, labels = sprintf("%.2f", zlim[1]),
                   cex = 0.7, col = "#444444", adj = c(0, 1))
    graphics::text(bx1, by0 - 0.015 * yr, labels = sprintf("%.2f", zlim[2]),
                   cex = 0.7, col = "#444444", adj = c(1, 1))
    graphics::text(mean(c(bx0, bx1)), by1 + 0.015 * yr,
                   labels = "net diversification", cex = 0.75,
                   col = "#111111", adj = c(0.5, 0))
  }

  invisible(x)
}
