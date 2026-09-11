# Whole-tree constant-rate birth-death fit.
#
# Thin wrapper around diversitree::make.bd() and find.mle(). Returns fitted
# lambda, mu, log-likelihood, and AIC, nothing else. No search, no
# candidates, this is the k = 0 model that every split model is compared
# against.

fit_whole_tree <- function(tree, sampling.f = 1, x.init = c(0.1, 0.03),
                            method = "subplex", lower = 0, upper = 10) {

  lik <- diversitree::make.bd(tree, sampling.f = sampling.f)
  fit <- diversitree::find.mle(lik, x.init, method = method,
                                lower = lower, upper = upper)

  pars <- stats::coef(fit)

  list(
    lambda = unname(pars["lambda"]),
    mu = unname(pars["mu"]),
    loglik = fit$lnLik,
    k = 2,
    AIC = 2 * 2 - 2 * fit$lnLik,
    fit = fit
  )
}
