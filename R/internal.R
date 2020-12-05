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

