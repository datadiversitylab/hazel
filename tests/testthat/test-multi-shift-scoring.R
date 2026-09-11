# Verify score_multi_shift_detection() against constructed cases with
# obvious known answers, using fake result and true_shifts objects so the
# scoring logic is checked in isolation from detection accuracy.

make_fake_tree <- function() {
  set.seed(1)
  tree <- ape::rtree(10)
  tree$tip.label <- paste0("t", 1:10)
  tree
}

test_that("perfect match scores as full recovery and exact match", {
  tree <- make_fake_tree()
  true_shifts <- list(list(shifted_clade_tips = c("t1", "t2", "t3")))
  mrca <- ape::getMRCA(tree, c("t1", "t2", "t3"))
  result <- list(splits = mrca)

  score <- score_multi_shift_detection(result, tree, true_shifts)
  expect_equal(score$n_recovered, 1)
  expect_true(score$exact_match)
})

test_that("two true shifts both detected score as full recovery, no false positives", {
  tree <- make_fake_tree()
  a <- c("t1", "t2"); b <- c("t7", "t8", "t9")
  true_shifts <- list(list(shifted_clade_tips = a),
                       list(shifted_clade_tips = b))
  result <- list(splits = c(ape::getMRCA(tree, a), ape::getMRCA(tree, b)))

  score <- score_multi_shift_detection(result, tree, true_shifts)
  expect_equal(score$n_recovered, 2)
  expect_equal(score$n_false_positive, 0)
})

test_that("one of two shifts detected scores as half recovery", {
  tree <- make_fake_tree()
  a <- c("t1", "t2"); b <- c("t7", "t8", "t9")
  true_shifts <- list(list(shifted_clade_tips = a),
                       list(shifted_clade_tips = b))
  result <- list(splits = ape::getMRCA(tree, a))

  score <- score_multi_shift_detection(result, tree, true_shifts)
  expect_equal(score$n_recovered, 1)
  expect_equal(score$fraction_recovered, 0.5)
  expect_false(score$all_recovered)
})

test_that("a spurious detection counts as a false positive without losing the true one", {
  tree <- make_fake_tree()
  a <- c("t1", "t2")
  true_shifts <- list(list(shifted_clade_tips = a))
  result <- list(splits = c(ape::getMRCA(tree, a),
                             ape::getMRCA(tree, c("t5", "t6"))))

  score <- score_multi_shift_detection(result, tree, true_shifts)
  expect_equal(score$n_recovered, 1)
  expect_equal(score$n_false_positive, 1)
})

test_that("null case returns NA fraction without error", {
  tree <- make_fake_tree()
  result <- list(splits = integer(0))
  score <- score_multi_shift_detection(result, tree, list())
  expect_equal(score$n_true, 0)
  expect_equal(score$n_detected, 0)
  expect_true(is.na(score$fraction_recovered))
})
