# Regression tests targeting specific properties verified by hand today,
# each with real, documented cost to re-diagnose from scratch if it ever
# silently breaks again.

testthat::test_that("fit_split_tree gives the same AIC regardless of node order in the nodes vector", {

  # diagnose_node77.R confirmed this by hand: {58, 77} and {77, 58} gave
  # the same AIC (both -225.29), the earlier "FALSE" result in that
  # script was a too-strict floating point tolerance, not a real order
  # effect. This makes that check permanent and uses a realistic
  # tolerance based on the optimizer's own convergence precision.

  set.seed(4)
  sim <- simulate_multi_shift_tree(
    n_background = 25,
    shifts = list(list(n_shifted = 10, lambda_ratio = 8),
                   list(n_shifted = 10, lambda_ratio = 6)),
    epsilon = 0.2
  )
  tree <- sim$tree
  candidates <- candidate_nodes(tree, min_clade_size = 4)
  testthat::skip_if(length(candidates) < 2, "not enough candidates for this check")

  two_nodes <- candidates[1:2]

  fit_forward <- fit_split_tree(tree, two_nodes, skip_reduction_check = TRUE)
  fit_reversed <- fit_split_tree(tree, rev(two_nodes), skip_reduction_check = TRUE)

  # subplex's own default reltol is on the order of 1e-4, two independent
  # optimizer runs should not be expected to agree past that
  testthat::expect_equal(fit_forward$AIC, fit_reversed$AIC, tolerance = 1e-2)
})

testthat::test_that("stepwise_clade_search's calibration mode returns NA gracefully with no candidates", {

  set.seed(1)
  tree <- simulate_null_tree(n = 10, rho = 1, epsilon = 0.2)

  step <- stepwise_clade_search(tree, candidates = integer(0),
                                 aic_threshold = NULL, verbose = FALSE)

  testthat::expect_true(is.na(step$best_node))
  testthat::expect_equal(step$aic_improvement, 0)
})

testthat::test_that("stepwise_clade_search's normal mode never accepts more than max_splits", {

  set.seed(1)
  tree <- simulate_null_tree(n = 40, rho = 1, epsilon = 0.2)
  candidates <- candidate_nodes(tree, min_clade_size = 4)
  testthat::skip_if(length(candidates) == 0, "no eligible candidates")

  result <- stepwise_clade_search(tree, candidates, aic_threshold = 0.01,
                                   max_splits = 2, verbose = FALSE)

  # threshold of 0.01 is deliberately permissive, so this checks the cap
  # itself is enforced, not whether the threshold is realistic
  testthat::expect_lte(length(result$splits), 2)
})

testthat::test_that("get_calibrated_threshold finds the exact matching row when one exists", {

  calib <- data.frame(n = c(20, 40, 60), rho = c(1, 1, 1),
                       epsilon = c(0.2, 0.2, 0.2), aic_threshold = c(10, 20, 30))

  result <- get_calibrated_threshold(calib, n = 40, rho = 1, epsilon = 0.2)
  testthat::expect_equal(result, 20)
})

testthat::test_that("candidate_nodes never proposes the root itself", {

  set.seed(1)
  tree <- simulate_null_tree(n = 20, rho = 1, epsilon = 0.2)
  root_node <- length(tree$tip.label) + 1

  candidates <- candidate_nodes(tree, min_clade_size = 2)
  testthat::expect_false(root_node %in% candidates)
})

testthat::test_that("candidate_nodes respects min_clade_size on both sides of the split", {

  set.seed(1)
  tree <- simulate_null_tree(n = 30, rho = 1, epsilon = 0.2)
  min_size <- 5

  candidates <- candidate_nodes(tree, min_clade_size = min_size)

  clade_sizes <- vapply(candidates, function(nd) {
    length(phangorn::Descendants(tree, nd, type = "tips")[[1]])
  }, integer(1))

  n_tips <- length(tree$tip.label)
  testthat::expect_true(all(clade_sizes >= min_size))
  testthat::expect_true(all(clade_sizes <= n_tips - min_size))
})
