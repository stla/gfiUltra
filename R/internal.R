#' @importFrom utils combn
#' @noRd
powerset <- function(ix, nms) {
  subsets <- unlist(lapply(seq_along(ix), function(i){
    combn(ix, i, simplify = FALSE)
  }), recursive = FALSE)
  names(subsets) <- vapply(
    subsets,
    function(subset) paste0(nms[subset], collapse = " - "),
    FUN.VALUE = character(1L)
  )
  subsets
}

#' @importFrom stats lm.fit
#' @noRd
RSS <- function(X, y, model){
  fit <- lm.fit(X[, model+1L, drop = FALSE], y)
  sum(fit[["residuals"]]^2)
}

modelsWithRSS <- function(X, y, models) {
  lapply(models, function(model){
    list(indices = model, rss = RSS(X, y, model))
  })
}

Rgamma_unnrmlzd <- function(n, p, modelWithRSS, gamma){
  d <- length(modelWithRSS[["indices"]]) + 1L
  logR <- lgamma((n-d)/2) - (n-d-1)/2 * log(pi * modelWithRSS[["rss"]]) -
    (d+1)/2 * log(n) - gamma * lchoose(p, d)
  exp(logR)
}

Rgammas <- function(n, p, models_with_RSS, gamma = 1){
  # n <- length(y)
  # p <- ncol(X)
  Rgammas_unnrmlzd <- vapply(
    models_with_RSS,
    function(model) Rgamma_unnrmlzd(n, p, model, gamma = 1),
    FUN.VALUE = numeric(1L)
  )
  Rgammas_unnrmlzd / sum(Rgammas_unnrmlzd)
}

