# Per-regime rate comparison: a horizontal lollipop chart of each
# regime's net diversification rate, with the background rate drawn as a
# reference line, so "which clades diversify faster or slower than
# background" is answerable at a glance. Points colored by the same rate
# ramp as the tree plot, keeping the visual language consistent.

plot_regime_rates <- function(x, ...) {

  if (!inherits(x, "rze_result")) stop("x must be an rze_result")

  rates <- x$rates
  if (nrow(rates) <= 1) {
    message("No shifts detected; nothing to compare against background.")
    return(invisible(NULL))
  }

  background_rate <- rates$net_diversification[rates$regime == "background"]

  # Order regimes by rate for a clean visual ranking, background last
  shift_rows <- rates[rates$regime != "background", ]
  shift_rows <- shift_rows[order(shift_rows$net_diversification), ]

  vals <- shift_rows$net_diversification
  labels <- sprintf("%s (n=%d)", shift_rows$regime, shift_rows$n_tips)

  zlim <- range(rates$net_diversification, na.rm = TRUE)
  pal <- rze_palette(256)
  point_cols <- rze_map_color(vals, palette = pal, zlim = zlim)

  old_par <- graphics::par(mar = c(4.5, 9, 3, 2), bg = "white")
  on.exit(graphics::par(old_par))

  y <- seq_along(vals)
  xlim <- range(c(vals, background_rate), na.rm = TRUE)
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.1

  graphics::plot(NA, xlim = xlim, ylim = c(0.5, length(vals) + 0.5),
                 yaxt = "n", xlab = "net diversification rate", ylab = "",
                 main = "Regime rates vs background",
                 col.main = "#111111", font.main = 1, cex.main = 1.3,
                 col.axis = "#444444", col.lab = "#222222")

  # Background reference line
  graphics::abline(v = background_rate, col = rze_neutral(), lwd = 2, lty = 2)
  graphics::text(background_rate, length(vals) + 0.4,
                 labels = sprintf("background (%.3f)", background_rate),
                 col = rze_neutral(), cex = 0.85, pos = 4, offset = 0.2)

  # Lollipop stems from background to each regime's rate
  graphics::segments(background_rate, y, vals, y, col = "#999999", lwd = 2)

  # Points
  graphics::points(vals, y, pch = 21, bg = point_cols, col = "#222222",
                   cex = 2.4, lwd = 1.5)

  # Regime labels on the y axis
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE,
                 col.axis = "#222222", cex.axis = 0.95)

  # Flag boundary_mu regimes with a small mark
  boundary_flag <- shift_rows$boundary_mu
  if (any(boundary_flag)) {
    graphics::points(vals[boundary_flag], y[boundary_flag], pch = 4,
                     col = "#222222", cex = 1.1, lwd = 2)
    graphics::mtext("x = extinction at boundary (rate less certain)",
                    side = 1, line = 3.3, cex = 0.75, col = "#666666", adj = 1)
  }

  invisible(NULL)
}
