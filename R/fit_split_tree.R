# Split-tree birth-death fit, one lambda and mu per partition.
#
# Thin wrapper around diversitree::make.bd.split() and find.mle(). nodes is
# a vector of node labels or ape node indices where the tree is split,
# split.t = Inf at every node, splitting at the base of the branch, the
# same convention MEDUSA uses.
#
# Every call runs one mandatory check before returning: fit the same tree
# with make.bd() alone, and confirm that evaluating the split likelihood
# at equal parameters across every partition reproduces the whole-tree
# likelihood. This is not a one-time test, it is a runtime assertion, run
# on every call, because a silent mismatch here (the wrong default, the
# wrong node indexing) is exactly the kind of thing that goes unnoticed
# for hours otherwise.

#' Fit a birth-death model with clade-specific regimes
#'
#' Fits one speciation and extinction pair per partition, splitting the tree
#' at the given nodes. Includes a runtime check that equal-rate partitions
#' reproduce the whole-tree likelihood.
#'
#' @param tree A \code{phylo} object.
#' @param nodes Node indices where the tree is split.
#' @param sampling.f Sampling fraction.
#' @param x.init Starting values per partition.
#' @param method Optimizer method.
#' @param lower,upper Parameter bounds.
#' @param skip_reduction_check Skip the equal-rate reduction check.
#' @return A list with per-partition rates, log-likelihood, AIC, and a
#'   boundary flag.
#' @export
fit_split_tree <- function(tree, nodes, sampling.f = 1,
                            x.init = c(0.1, 0.03), method = "subplex",
                            lower = 0, upper = 10,
                            skip_reduction_check = FALSE) {

  n_partitions <- length(nodes) + 1

  lik_split <- diversitree::make.bd.split(
    tree, nodes,
    split.t = rep(Inf, length(nodes)),
    sampling.f = sampling.f
  )

  if (!skip_reduction_check) {
    lik_whole <- diversitree::make.bd(tree, sampling.f = sampling.f)
    equal_pars <- rep(x.init, n_partitions)
    names(equal_pars) <- diversitree::argnames(lik_split)

    split_at_equal <- lik_split(equal_pars)
    whole_at_same <- lik_whole(x.init)

    if (abs(split_at_equal - whole_at_same) > 1e-3) {
      stop(sprintf(
        paste("Split model does not reduce to the whole-tree model at",
              "equal parameters (split = %.5f, whole = %.5f). Something",
              "is wrong with the node specification or the engine before",
              "this fit can be trusted."),
        split_at_equal, whole_at_same
      ))
    }
  }

  x.init.full <- rep(x.init, n_partitions)
  names(x.init.full) <- diversitree::argnames(lik_split)

  fit <- diversitree::find.mle(lik_split, x.init.full, method = method,
                                lower = lower, upper = upper)
  pars <- stats::coef(fit)

  lambda <- pars[grepl("^lambda", names(pars))]
  mu <- pars[grepl("^mu", names(pars))]

  k <- length(pars)

  list(
    lambda = unname(lambda),
    mu = unname(mu),
    boundary_mu = any(mu < 1e-6),
    loglik = fit$lnLik,
    k = k,
    AIC = 2 * k - 2 * fit$lnLik,
    nodes = nodes,
    fit = fit
  )
}
