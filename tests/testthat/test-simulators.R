# Unit tests for simulate_null_tree, simulate_clade_shift_tree, and
# simulate_multi_shift_tree. Several of these tests exist specifically
# because they would have caught bugs found by hand today: a graft that
# silently broke ultrametricity (caught only by manually inspecting
# root-to-tip distances), and a tip-label collision between the
# background and donor trees (caught only by an error message during an
# unrelated test run). Both should now fail loudly and immediately if
# they ever recur, not require another afternoon of manual diagnosis.

testthat::test_that("simulate_null_tree produces a valid ultrametric tree of the requested size", {

  set.seed(1)
  tree <- simulate_null_tree(n = 25, rho = 1, epsilon = 0.2)

  testthat::expect_s3_class(tree, "phylo")
  testthat::expect_equal(length(tree$tip.label), 25)
  testthat::expect_true(ape::is.ultrametric(tree))
})

testthat::test_that("simulate_clade_shift_tree produces an ultrametric tree", {

  set.seed(1)
  sim <- simulate_clade_shift_tree(n_background = 20, n_shifted = 12,
                                    lambda_ratio = 5, epsilon = 0.2)

  # This exact check failed silently three times today before the
  # underlying edge-splitting bug was found and fixed. It stays a
  # dedicated, explicit test rather than folding into a larger one, so a
  # future regression here is immediately unambiguous about what broke.
  testthat::expect_true(ape::is.ultrametric(sim$tree))
})

testthat::test_that("simulate_clade_shift_tree's shifted clade tips are monophyletic in the grafted tree", {

  set.seed(1)
  sim <- simulate_clade_shift_tree(n_background = 20, n_shifted = 12,
                                    lambda_ratio = 5, epsilon = 0.2)

  mrca_node <- ape::getMRCA(sim$tree, sim$shifted_clade_tips)
  descendant_tips <- sim$tree$tip.label[
    phangorn::Descendants(sim$tree, mrca_node, type = "tips")[[1]]
  ]

  testthat::expect_setequal(descendant_tips, sim$shifted_clade_tips)
})

testthat::test_that("simulate_clade_shift_tree never produces duplicated tip labels", {

  # TreeSim assigns generic labels (t1, t2, ...) independently to each
  # simulated tree, background and donor clades collided on names before
  # this was fixed with an explicit relabeling step. Run several seeds,
  # since a collision depends on which specific labels each simulation
  # happens to draw.
  for (s in 1:10) {
    set.seed(s)
    sim <- simulate_clade_shift_tree(n_background = 15, n_shifted = 10,
                                      lambda_ratio = 4, epsilon = 0.2)
    testthat::expect_false(any(duplicated(sim$tree$tip.label)),
                            info = sprintf("seed %d produced duplicate labels", s))
  }
})

testthat::test_that("simulate_multi_shift_tree produces an ultrametric tree with multiple grafts", {

  set.seed(1)
  sim <- simulate_multi_shift_tree(
    n_background = 25,
    shifts = list(list(n_shifted = 10, lambda_ratio = 8),
                   list(n_shifted = 10, lambda_ratio = 6)),
    epsilon = 0.2
  )

  testthat::expect_true(ape::is.ultrametric(sim$tree))
  testthat::expect_false(any(duplicated(sim$tree$tip.label)))
  testthat::expect_length(sim$shifts, 2)
})

testthat::test_that("simulate_multi_shift_tree's shifted clades do not overlap with each other", {

  set.seed(1)
  sim <- simulate_multi_shift_tree(
    n_background = 25,
    shifts = list(list(n_shifted = 10, lambda_ratio = 8),
                   list(n_shifted = 10, lambda_ratio = 6)),
    epsilon = 0.2
  )

  overlap <- intersect(sim$shifts[[1]]$shifted_clade_tips,
                        sim$shifts[[2]]$shifted_clade_tips)
  testthat::expect_length(overlap, 0)
})
