# Tests for resolve_calibration() and the table validation/save helpers.

make_valid_table <- function() {
  data.frame(n = c(20, 50, 100), rho = 1, epsilon = 0.2,
             aic_threshold = c(10, 15, 20))
}

test_that("a data frame passed directly is returned as-is", {
  tab <- make_valid_table()
  got <- resolve_calibration(tab, verbose = FALSE)
  expect_equal(got, tab)
})

test_that("a valid table passes validation, an invalid one is rejected", {
  expect_true(validate_calibration_table(make_valid_table()))

  bad <- make_valid_table()
  bad$aic_threshold <- NULL
  expect_error(validate_calibration_table(bad), "missing required column")

  all_na <- make_valid_table()
  all_na$aic_threshold <- NA_real_
  expect_error(validate_calibration_table(all_na), "no usable thresholds")
})

test_that("a table saved to .rds is loaded back correctly", {
  tab <- make_valid_table()
  path <- tempfile(fileext = ".rds")
  saveRDS(tab, path)

  got <- resolve_calibration(path, verbose = FALSE)
  expect_equal(got, tab)
  unlink(path)
})

test_that("a table saved to .rda is loaded back correctly", {
  calibration_default <- make_valid_table()
  path <- tempfile(fileext = ".rda")
  save(calibration_default, file = path)

  got <- resolve_calibration(path, verbose = FALSE)
  expect_equal(got, calibration_default)
  unlink(path)
})

test_that("a nonexistent file path errors clearly", {
  expect_error(resolve_calibration("/no/such/file.rds", verbose = FALSE),
               "not found")
})

test_that("NULL returns NULL when no default is shipped, so caller can fall back", {
  # In the test environment the package data may or may not be installed;
  # this only checks that NULL input does not error and yields NULL or a
  # valid table, never a crash
  got <- tryCatch(resolve_calibration(NULL, verbose = FALSE),
                  error = function(e) "errored")
  expect_true(is.null(got) || is.data.frame(got))
})

test_that("save_as_default_calibration writes an .rda into data/", {
  tab <- make_valid_table()
  root <- tempfile()
  dir.create(root)

  path <- save_as_default_calibration(tab, package_root = root)
  expect_true(file.exists(path))
  expect_match(path, "calibration_default\\.rda$")

  unlink(root, recursive = TRUE)
})
