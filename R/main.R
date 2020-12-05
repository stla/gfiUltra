#' Title
#' @description xxx
#'
#' @param formula xxx
#' @param data xxx
#' @param gamma xxx
#' @param ... named arguments passed to \code{\link[SIS:SIS]{SIS}}
#'
#' @return xxx
#'
#' @export
#' @importFrom SIS SIS
#' @importFrom lazyeval f_eval_lhs
#' @importFrom stats model.matrix rchisq
#' @importFrom mvtnorm rmvnorm
#'
#' @examples xxx
gfiUltra <- function(formula, data, nsims = 5L, gamma = 1, ...){
  y <- f_eval_lhs(formula, data = data)
  X <- model.matrix(formula, data = data)
  n <- length(y)
  p <- ncol(X)
  Xm1 <- X[, -1L, drop = FALSE]
  sis <- suppressMessages(SIS(Xm1, y = y, family = "gaussian", ...))
  models <- powerset(sis[["ix"]], colnames(Xm1))
  models_with_FIT <- modelsWithFIT(X = X, y = y, models = models)
  modelsProbs <-
    Rgammas(n = n, p = p, models_with_FIT = models_with_FIT, gamma = gamma)
  #
  # simulations
  selected <- colnames(X)[c(1L, 1L + sis[["ix"]])]
  Sims <- matrix(0, nrow = nsims, ncol = length(selected) + 1L)
  colnames(Sims) <- c(selected, "sigma")
  nmodels <- length(models)
  for(i in 1L:nsims){
    M <- models_with_FIT[[sample.int(nmodels, 1L, prob = modelsProbs)]]
    Mvariables <- c("(Intercept)", names(M[["indices"]])) # donc on s'en fout des indices
    Mdim <- length(Mvariables)
    Mrss <- M[["fit"]][["rss"]]
    sigma2 <- Mrss / rchisq(1L, n - Mdim)
    Mcoeffs <- M[["fit"]][["coeffs"]]
    Sigma <- sigma2*M[["fit"]][["XtXinv"]]
    beta <- rmvnorm(1L, mean = Mcoeffs, sigma = Sigma, checkSymmetry = FALSE)
    Sims[i, c(Mvariables, "sigma")] <- c(beta, sqrt(sigma2))
  }
  Sims
}
