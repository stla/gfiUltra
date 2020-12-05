#' @importFrom utils combn
#' @noRd
powerset <- function(ix, nms) {
  subsets <- unlist(lapply(seq_along(ix), function(i){
    combn(ix, i, simplify = FALSE)
  }), recursive = FALSE)
  subsets <- lapply(subsets, function(subset){
    `names<-`(subset, nms[subset])
  })
  names(subsets) <- vapply(
    subsets,
    function(subset) paste0(names(subset), collapse = " - "),
    FUN.VALUE = character(1L)
  )
  subsets
}

#' @importFrom stats lm.fit
#' @noRd
lmFit <- function(X, y){
  fit <- lm.fit(X, y)
  list(
    rss = sum(fit[["residuals"]]^2L),
    coeffs = fit[["coefficients"]],
    XtXinv = tcrossprod(qr.solve(fit[["qr"]], diag(length(y))))
  )
}

FITmodel <- function(X, y, model){
  lmFit(X[, c(1L, 1L + model), drop = FALSE], y)
}

modelsWithFIT <- function(X, y, models) {
  lapply(models, function(model){
    list(indices = model, fit = FITmodel(X, y, model))
  })
}

Rgamma_unnrmlzd <- function(n, p, modelWithFIT, gamma){
  d <- length(modelWithFIT[["indices"]]) + 1L
  logR <- lgamma((n-d)/2) - (n-d-1)/2 * log(pi * modelWithFIT[["fit"]][["rss"]]) -
    (d+1)/2 * log(n) - gamma * lchoose(p, d)
  exp(logR)
}

Rgammas <- function(n, p, models_with_FIT, gamma){
  # n <- length(y)
  # p <- ncol(X)
  Rgammas_unnrmlzd <- vapply(
    models_with_FIT,
    function(model) Rgamma_unnrmlzd(n, p, model, gamma),
    FUN.VALUE = numeric(1L)
  )
  Rgammas_unnrmlzd / sum(Rgammas_unnrmlzd)
}

