# rze

Calibrated detection of clade-localized diversification shifts on time-calibrated phylogenies of extant taxa.

## What it does

rze finds where on a phylogeny net diversification rate shifts, which clades diversify faster or slower than the rest of the tree, and reports the net diversification rate for each detected regime. It uses a stepwise model-selection search over candidate clades, built on the birth-death split-model likelihood from `diversitree` (the same likelihood MEDUSA uses), with one important addition: the AIC threshold for accepting a shift is calibrated against simulated constant-rate trees, so the false-positive rate is controlled rather than left to an uncalibrated default.

## Why net diversification, and not separate speciation and extinction

For an extant-only tree, speciation and extinction are not separately identifiable in general (Louca & Pennell 2020). Net diversification is far better constrained, and it is what rze reports as its primary output. The underlying model still fits speciation and extinction internally, and rze flags any regime where the extinction estimate hit its boundary, so you always know when the separation is unsupported by the data.

## Quick start

```r
library(rze)

# your tree: an ultrametric phylo object of extant taxa
result <- rze(tree, rho = 0.8)   # rho = fraction of species sampled

print(result)     # detected shifts and per-regime net diversification
plot(result)      # branch-colored tree, shifts marked
plot_regime_rates(result)   # each regime's rate vs background
```

Supplying a calibration table (recommended for speed and rigor):

```r
calib <- calibrate_rze(n_grid = length(tree$tip.label),
                       rho_grid = 0.8, epsilon_grid = 0.2,
                       n_replicates = 100)
result <- rze(tree, calibration = calib, rho = 0.8)
```

## Visuals

rze ships bold, poster-ready base R plots: a branch-colored phylogeny where every branch is colored by its regime's net diversification rate, and a per-regime comparison showing which clades stand out against background. The visual language is consistent across every plot.

## Honest scope

- Detection power depends strongly on tree size and shift magnitude. Small trees (roughly under 30 tips) have low power even for large shifts; recovery is reliable for large trees. See the `rze-benchmarks` repository for full power curves with confidence intervals.
- Sampling fraction is currently global. Per-clade sampling is a documented future extension; the empirical harness in `rze-benchmarks` shows how to handle sampling uncertainty robustly in the meantime.
- rze reports net diversification, not separately identified speciation and extinction.

## License

MIT
