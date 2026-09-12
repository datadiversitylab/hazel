# hazel

Calibrated detection of clade-localized diversification shifts on time-calibrated phylogenies of extant taxa.

## What it does

hazel finds where net diversification rate shifts on a phylogeny, which clades diversify faster or slower than the rest of the tree, and reports the net diversification rate for each detected regime. It uses a stepwise search over candidate clades with an AIC threshold calibrated against simulated constant-rate trees, so the false-positive rate is controlled rather than left to an uncalibrated default. A pre-computed calibration table ships with the package, so most analyses need no calibration step.

## Install and build

```r
library(devtools)
install_github("datadiversitylab/hazel")
```

## Quick start

```r
library(hazel)

result <- hazel(my_tree, rho = 0.8)   # rho = fraction of species sampled

result                # detected shifts and per-regime net diversification
plot(result)          # branch-colored tree, shifts marked
plot_regime_rates(result)   # each regime's rate vs background
```

The shipped calibration table is used automatically. You only run `calibrate_hazel()` yourself if your tree falls outside the shipped grid, or you want a threshold tuned to your exact settings.

## Vignettes

- **Getting started with hazel** — the common case, using the shipped table.
- **Calibration: what it does and when you need it** — what the calibration table is and how to build your own.
- **Two ways to run hazel** — the shipped table and your own calibration, side by side.

Build them with:

```r
devtools::build_vignettes()
```

## Scope, stated plainly

- Detection power depends on tree size and shift magnitude. Small trees (under ~30 tips) have low power even for large shifts; recovery is reliable for large trees.
- Sampling fraction is currently global. Uneven sampling is handled by running across plausible values and keeping shifts robust across them.
- hazel reports net diversification, not separately identified speciation and extinction, because only net diversification is reliably identifiable from an extant-only tree.

## License

MIT
