# Resolving which calibration table to use.
#
# rze needs an AIC-threshold calibration table to run. There are three
# ways one can be supplied, in precedence order:
#
#   1. An explicit table passed by the user (a data frame, or a path to a
#      saved .rds/.rda file) -- always wins.
#   2. The default table shipped with the package (data/calibration_default.rda),
#      loaded automatically if present.
#   3. On-the-fly calibration for the specific tree, the slow fallback,
#      used only when neither of the above is available.
#
# This function returns a calibration data frame given whatever the user
# passed (which may be NULL, a data frame, or a file path), falling back
# through the order above.

resolve_calibration <- function(calibration = NULL, verbose = TRUE) {

  # Case 1a: user passed a data frame directly
  if (is.data.frame(calibration)) {
    validate_calibration_table(calibration)
    return(calibration)
  }

  # Case 1b: user passed a path to a saved table
  if (is.character(calibration) && length(calibration) == 1) {
    if (!file.exists(calibration)) {
      stop("Calibration file not found: ", calibration)
    }
    tab <- read_calibration_file(calibration)
    validate_calibration_table(tab)
    if (verbose) message("Loaded calibration table from ", calibration)
    return(tab)
  }

  # Case 2: shipped default, if the package bundles one
  default <- tryCatch({
    e <- new.env()
    utils::data("calibration_default", package = "rze", envir = e)
    get("calibration_default", envir = e)
  }, error = function(err) NULL, warning = function(w) NULL)

  if (!is.null(default)) {
    validate_calibration_table(default)
    if (verbose) message("Using the calibration table shipped with rze.")
    return(default)
  }

  # Case 3: signal that no table is available, caller decides whether to
  # calibrate on the fly
  NULL
}

# Read a calibration table from either an .rds or an .rda file
read_calibration_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    return(readRDS(path))
  }
  if (ext == "rda" || ext == "rdata") {
    e <- new.env()
    loaded <- load(path, envir = e)
    # take the first data frame found in the loaded objects
    for (nm in loaded) {
      obj <- get(nm, envir = e)
      if (is.data.frame(obj)) return(obj)
    }
    stop("No data frame found in ", path)
  }
  stop("Calibration file must be .rds or .rda, got: ", ext)
}

# Confirm a calibration table has the columns the lookup needs
validate_calibration_table <- function(tab) {
  required <- c("n", "rho", "epsilon", "aic_threshold")
  missing <- setdiff(required, names(tab))
  if (length(missing) > 0) {
    stop("Calibration table is missing required column(s): ",
         paste(missing, collapse = ", "),
         ". A valid table has columns: ", paste(required, collapse = ", "), ".")
  }
  if (all(is.na(tab$aic_threshold))) {
    stop("Calibration table has no usable thresholds (all NA).")
  }
  invisible(TRUE)
}

# Save a generated calibration table as the package's shipped default.
# Run this once, on real hardware, after generating a production-grade
# table with calibrate_rze() across a full grid. Writes into the package
# source tree's data/ directory so it ships on the next install/build.
save_as_default_calibration <- function(calibration_table,
                                         package_root = ".") {
  validate_calibration_table(calibration_table)
  data_dir <- file.path(package_root, "data")
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
  calibration_default <- calibration_table
  save(calibration_default,
       file = file.path(data_dir, "calibration_default.rda"),
       compress = "xz")
  message("Saved shipped calibration to ",
          file.path(data_dir, "calibration_default.rda"),
          " (", nrow(calibration_table), " cells).")
  invisible(file.path(data_dir, "calibration_default.rda"))
}
