# Shared visual identity for rze plots. Every plotting function draws from
# this so the package has one coherent, bold, talk-and-poster-ready look
# rather than ad hoc styling per plot. Base R graphics throughout.
#
# The palette is a single-hue, high-saturation sequential ramp from a
# deep cool low to a hot high, chosen to read clearly from across a room
# and to reproduce well in projection. It is perceptually ordered (dark
# low to bright high), avoiding the rainbow ramp that misleads the eye and
# dates a figure.

rze_palette <- function(n = 256) {
  # Deep indigo -> violet -> magenta -> hot orange -> bright yellow, a
  # saturated sequential ramp with strong luminance ordering
  anchors <- c("#1a1040", "#5b1d8c", "#a11a9c", "#e03b8b",
               "#f96e3d", "#ffb020", "#ffe94d")
  grDevices::colorRampPalette(anchors)(n)
}

# Map a numeric vector to palette colors over a given range
rze_map_color <- function(values, palette = rze_palette(256),
                          zlim = range(values, na.rm = TRUE)) {
  if (diff(zlim) == 0) zlim <- zlim + c(-1, 1) * 0.5
  idx <- round((values - zlim[1]) / diff(zlim) * (length(palette) - 1)) + 1
  idx[idx < 1] <- 1
  idx[idx > length(palette)] <- length(palette)
  cols <- palette[idx]
  cols[is.na(values)] <- "#cccccc"
  cols
}

# Consistent graphical parameters, applied inside each plot function via
# on.exit-restored par(). Dark-friendly, generous, clean.
rze_par <- function() {
  list(
    bg = "white",
    fg = "#222222",
    font.main = 1,
    col.main = "#111111",
    cex.main = 1.3,
    col.axis = "#444444",
    col.lab = "#222222",
    family = ""  # left to the device default, avoids font-not-found issues
  )
}

# The accent color for marking detected shifts, chosen to sit off the
# rate ramp so shift markers never blend into a branch color
rze_shift_color <- function() "#00d0c0"

# A neutral color for background/context elements
rze_neutral <- function() "#333333"

# Draw a horizontal color bar legend for the rate scale. Standalone so it
# can be placed in its own panel or beside a tree.
rze_colorbar <- function(zlim, palette = rze_palette(256),
                         label = "net diversification",
                         horiz = TRUE) {

  old <- graphics::par(mar = c(3, 1, 2, 1))
  on.exit(graphics::par(old))

  n <- length(palette)
  if (horiz) {
    graphics::plot.new()
    graphics::plot.window(xlim = zlim, ylim = c(0, 1))
    xs <- seq(zlim[1], zlim[2], length.out = n + 1)
    graphics::rect(xs[-(n + 1)], 0, xs[-1], 1, col = palette, border = NA)
    graphics::rect(zlim[1], 0, zlim[2], 1, border = rze_neutral(), lwd = 1.5)
    graphics::axis(1, col = rze_neutral(), col.axis = "#444444", cex.axis = 0.9)
    graphics::title(main = label, cex.main = 1.0, col.main = "#111111", font.main = 1)
  }
  invisible(NULL)
}
