# Internal nodes eligible as split points, the candidate set for
# stepwise_clade_search(). A node is eligible if the clade it subtends has
# at least min_clade_size tips, so a search can't propose partitioning off
# a handful of tips with almost no data behind the estimate.

candidate_nodes <- function(tree, min_clade_size = 5) {

  n_tips <- length(tree$tip.label)
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)

  clade_sizes <- vapply(internal_nodes, function(nd) {
    length(phangorn::Descendants(tree, nd, type = "tips")[[1]])
  }, integer(1))

  # Exclude the root itself, it is already the whole-tree partition
  root_node <- n_tips + 1
  eligible <- internal_nodes[internal_nodes != root_node &
                                clade_sizes[internal_nodes - n_tips] >= min_clade_size &
                                clade_sizes[internal_nodes - n_tips] <= n_tips - min_clade_size]

  eligible
}
