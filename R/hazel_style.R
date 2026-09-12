# The hazel visual identity: an earthy palette and a small set of shared
# drawing helpers, so every hazel figure has one calm, recognizable look.
# Colors run from cool (slow) to warm (fast): teal and green through gold
# and orange to deep rust. Both ends carry full color, no black, no white.

#' The hazel color palette
#'
#' An earthy sequential ramp used across hazel plots, running from cool
#' teal-green tones (low rates) through gold and orange to deep rust (high
#' rates). Both ends are richly saturated, so the fastest and slowest
#' clades read as clearly as everything in between. Suitable for shading
#' branches or filling a rate legend.
#'
#' @param n Number of colors to return.
#' @return A character vector of hex colors.
#' @export
hazel_palette <- function(n = 256) {
  anchors <- c("#2f5d62", "#5c8a5a", "#c9b458", "#e79a3c",
               "#d76b2c", "#b83f22", "#7a2418")
  grDevices::colorRampPalette(anchors)(n)
}

# Distinct colors for separate regimes (clades), drawn from the same earthy
# family so they sit together without clashing.
hazel_regime_colors <- function() {
  c("#3a7a5e", "#d76b2c", "#c9a13e", "#b83f22",
    "#2f5d62", "#8a6f4e", "#5c8a5a", "#7a2418")
}

# Map numeric values onto the palette over a given range.
hazel_map_color <- function(values, palette = hazel_palette(256),
                            zlim = range(values, na.rm = TRUE)) {
  if (diff(zlim) == 0) zlim <- zlim + c(-1, 1) * 0.5
  idx <- round((values - zlim[1]) / diff(zlim) * (length(palette) - 1)) + 1
  idx[idx < 1] <- 1
  idx[idx > length(palette)] <- length(palette)
  cols <- palette[idx]
  cols[is.na(values)] <- "#cabfa8"
  cols
}

# The color used to mark detected shifts: a soft cream ring that reads
# clearly against every band of the earthy ramp.
hazel_mark <- function() "#f4ecdd"
hazel_mark_border <- function() "#2b1a12"

# A quiet ink color for axes, labels, and fine lines.
hazel_ink <- function() "#3d3d3d"

# Draw a compact per-regime key: a small color swatch, the regime name, and
# its net diversification rate, one row per regime, in the lower left. No
# plot title, in keeping with the hazel look. Regimes whose extinction hit
# the boundary are marked so the reader knows that rate is less certain.
hazel_regime_key <- function(regime_names, rates, colors, boundary = NULL) {
  usr <- graphics::par("usr")
  xr <- usr[2] - usr[1]; yr <- usr[4] - usr[3]
  
  n <- length(regime_names)
  # anchor rows in the lower-left, stacking upward
  x0 <- usr[1] + 0.03 * xr
  sw <- 0.022 * xr                       # swatch width
  row_h <- 0.045 * yr                    # row spacing
  y_base <- usr[3] + 0.03 * yr
  
  for (i in seq_len(n)) {
    yy <- y_base + (n - i) * row_h
    graphics::rect(x0, yy, x0 + sw, yy + 0.028 * yr,
                   col = colors[i], border = hazel_ink(), lwd = 0.6)
    label <- sprintf("%s   %.3f", regime_names[i], rates[i])
    if (!is.null(boundary) && isTRUE(boundary[i])) label <- paste0(label, " *")
    graphics::text(x0 + sw + 0.012 * xr, yy + 0.014 * yr, labels = label,
                   cex = 0.7, col = hazel_ink(), adj = c(0, 0.5))
  }
  
  graphics::text(x0, y_base + n * row_h + 0.008 * yr,
                 labels = "net diversification", cex = 0.72,
                 col = hazel_ink(), adj = c(0, 0))
  if (!is.null(boundary) && any(boundary)) {
    graphics::text(x0, y_base - 0.028 * yr,
                   labels = "* extinction at boundary; rate less certain",
                   cex = 0.62, col = "#8a8175", adj = c(0, 1))
  }
  invisible(NULL)
}