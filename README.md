# rze

Calibrated detection of clade-localized diversification shifts on time-calibrated phylogenies of extant taxa.

## What it does

rze finds where on a phylogeny net diversification rate shifts, which clades diversify faster or slower than the rest of the tree, and reports the net diversification rate for each detected regime. It uses a stepwise model-selection search over candidate clades, built on the birth-death split-model likelihood from `diversitree` (the same likelihood MEDUSA uses), with one important addition: the AIC threshold for accepting a shift is calibrated against simulated constant-rate trees, so the false-positive rate is controlled rather than left to an uncalibrated default.

## Why net diversification, and not separate speciation and extinction

For an extant-only tree, speciation and extinction are not separately identifiable in general (Louca & Pennell 2020). Net diversification is far better constrained, and it is what rze reports as its primary output. The underlying model still fits speciation and extinction internally, and rze flags any regime where the extinction estimate hit its boundary, so you always know when the separation is unsupported by the data.

## Quick start

Install rze from GitHub

```
library(devtools)
install_github("datadiversitylab/rze")
```

Below is a brief and quick guide to rze.

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

## Calibration tables

rze needs an AIC-threshold calibration table. It resolves one in this order:

1. **A table you pass in** — either a data frame from `calibrate_rze()`, or a path to a saved `.rds`/`.rda` file. Always wins.
2. **The table shipped with the package**, loaded automatically if present. Generate it once with `scripts/generate_shipped_calibration.R` in the `rze-benchmarks` repository, which runs a full-grid, high-replicate calibration and saves it into the package via `save_as_default_calibration()`.
3. **On-the-fly calibration** for your specific tree — the slow fallback, used only when neither of the above is available.

To use a pre-existing table you generated:

```r
result <- rze(tree, calibration = "my_calibration.rds", rho = 0.8)
```

To make a generated table the package's shipped default:

```r
save_as_default_calibration(my_table, package_root = "path/to/rze")
# then rebuild/reinstall the package
```

## Visuals

rze ships bold, poster-ready base R plots: a branch-colored phylogeny where every branch is colored by its regime's net diversification rate, and a per-regime comparison showing which clades stand out against background. The visual language is consistent across every plot.

## Honest scope

- Detection power depends strongly on tree size and shift magnitude. Small trees (roughly under 30 tips) have low power even for large shifts; recovery is reliable for large trees. See the `rze-benchmarks` repository for full power curves with confidence intervals.
- Sampling fraction is currently global. Per-clade sampling is a documented future extension; the empirical harness in `rze-benchmarks` shows how to handle sampling uncertainty robustly in the meantime.
- rze reports net diversification, not separately identified speciation and extinction.

## License

MIT
