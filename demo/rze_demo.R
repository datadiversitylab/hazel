# End-to-end demonstration of rze on a simulated tree with known shifts.
# Runs detection, prints the result, and produces every package plot.
# Also serves as the source for the README figures.

library(ape)
library(diversitree)
library(phangorn)
library(TreeSim)
library(phytools)

# Source all package functions (for running before formal install)
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

set.seed(42)

# A tree with three real, clearly-different clade shifts, sized large
# enough to sit in the reliable-detection region established by the
# benchmarks (recovery is high for multi-shift trees around this size).
sim <- simulate_multi_shift_tree(
  n_background = 150,
  shifts = list(
    list(n_shifted = 30, lambda_ratio = 8),
    list(n_shifted = 25, lambda_ratio = 5),
    list(n_shifted = 20, lambda_ratio = 3)
  ),
  epsilon = 0.2, rho = 1
)
tree <- sim$tree
cat(sprintf("Simulated tree: %d tips, 3 true shifts\n", length(tree$tip.label)))

# Calibrate for this tree size (in practice, use a shipped table or a
# larger calibration; 50 replicates here for a runnable demo)
cat("Calibrating...\n")
calib <- calibrate_rze(
  n_grid = length(tree$tip.label),
  rho_grid = 1, epsilon_grid = 0.2,
  n_replicates = 50, min_clade_size = 5, verbose = FALSE
)

# Run rze
result <- rze(tree, calibration = calib, rho = 1, epsilon = 0.2,
              min_clade_size = 5, verbose = TRUE)

print(result)

# Plot 1: the hero image, branch-colored tree
png("demo/rze_tree.png", width = 1000, height = 1100, res = 130)
plot(result)
dev.off()

# Plot 2: per-regime rate comparison
png("demo/rze_regime_rates.png", width = 1000, height = 700, res = 130)
plot_regime_rates(result)
dev.off()

cat("\nFigures written to demo/. Open rze_tree.png and rze_regime_rates.png.\n")
