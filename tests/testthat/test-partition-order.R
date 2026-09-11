# Pin down the partition-ordering convention that regime_rates() relies
# on. make.bd.split(tree, nodes) must order fitted parameters so that
# partition 1 is the background (root) regime and partition i+1
# corresponds to nodes[i]. regime_rates() assumes exactly this. If
# diversitree ever changes it, this test fails loudly rather than letting
# rates be silently attached to the wrong clades.
#
# The check: construct a tree with one clade that has a genuinely, obviously
# different rate (a dense burst), split at that clade's node, and confirm
# the partition whose fitted net diversification is highest is the one
# corresponding to that node, not the background partition.

test_that("make.bd.split partition order maps partition i+1 to nodes[i]", {

  set.seed(1)
  # Background tree with a fast clade grafted on, known node
  sim <- simulate_clade_shift_tree(n_background = 40, n_shifted = 25,
                                    lambda_ratio = 10, epsilon = 0.1)
  tree <- sim$tree
  fast_node <- ape::getMRCA(tree, sim$shifted_clade_tips)

  fit <- fit_split_tree(tree, fast_node, skip_reduction_check = TRUE)

  # fit$lambda[1] is background, fit$lambda[2] should be the fast clade.
  # Net diversification of partition 2 should exceed partition 1.
  net_div <- fit$lambda - fit$mu
  expect_length(net_div, 2)
  expect_gt(net_div[2], net_div[1])
})

test_that("fit_split_tree stores the node vector it was given, for downstream mapping", {
  set.seed(2)
  tree <- simulate_null_tree(n = 40, rho = 1, epsilon = 0.2)
  cands <- candidate_nodes(tree, min_clade_size = 5)
  skip_if(length(cands) < 2, "not enough candidates")

  chosen <- sort(cands[1:2])
  fit <- fit_split_tree(tree, chosen, skip_reduction_check = TRUE)
  expect_equal(fit$nodes, chosen)
})
