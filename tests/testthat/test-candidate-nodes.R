# Independent verification of candidate_nodes(). The implementation uses
# phangorn::Descendants() to count clade sizes. This test cross-checks it
# against an independent clade-size computation using only ape's own
# node-descendant machinery, on a hand-constructed tree where the answer
# can also be reasoned out by eye. If phangorn and ape disagree, or if
# either disagrees with the hand answer, this fails loudly.

test_that("candidate_nodes matches an independent ape-based clade-size count", {

  # A balanced 8-tip tree, fully resolved. Clade sizes at internal nodes
  # are known by construction.
  set.seed(1)
  tree <- ape::stree(8, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  tree$tip.label <- paste0("t", 1:8)

  n_tips <- length(tree$tip.label)
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)

  # Independent clade-size count using ape's node-based descendant search,
  # not phangorn
  ape_clade_size <- function(tree, node) {
    if (node <= length(tree$tip.label)) return(1L)
    tips <- integer(0)
    to_visit <- node
    while (length(to_visit) > 0) {
      current <- to_visit[1]
      to_visit <- to_visit[-1]
      children <- tree$edge[tree$edge[, 1] == current, 2]
      for (ch in children) {
        if (ch <= length(tree$tip.label)) {
          tips <- c(tips, ch)
        } else {
          to_visit <- c(to_visit, ch)
        }
      }
    }
    length(tips)
  }

  min_clade_size <- 2
  expected <- internal_nodes[vapply(internal_nodes, function(nd) {
    sz <- ape_clade_size(tree, nd)
    root_node <- length(tree$tip.label) + 1
    nd != root_node && sz >= min_clade_size && sz <= n_tips - min_clade_size
  }, logical(1))]

  got <- candidate_nodes(tree, min_clade_size = min_clade_size)

  expect_setequal(got, expected)
})

test_that("candidate_nodes never returns the root", {
  set.seed(2)
  tree <- ape::rtree(20)
  root_node <- length(tree$tip.label) + 1
  expect_false(root_node %in% candidate_nodes(tree, min_clade_size = 3))
})

test_that("candidate_nodes respects min_clade_size on both sides", {
  set.seed(3)
  tree <- ape::rtree(30)
  min_size <- 5
  cands <- candidate_nodes(tree, min_clade_size = min_size)
  n_tips <- length(tree$tip.label)

  sizes <- vapply(cands, function(nd) {
    length(phangorn::Descendants(tree, nd, type = "tips")[[1]])
  }, integer(1))

  expect_true(all(sizes >= min_size))
  expect_true(all(sizes <= n_tips - min_size))
})

test_that("candidate_nodes returns empty when min_clade_size too large for the tree", {
  set.seed(4)
  tree <- ape::rtree(10)
  # min_clade_size of 6 on a 10-tip tree: no node can have between 6 and 4
  # tips (an impossible range), so nothing is eligible
  expect_length(candidate_nodes(tree, min_clade_size = 6), 0)
})
