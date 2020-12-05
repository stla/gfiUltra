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
#' @importFrom stats model.matrix rchisq rnorm
#'
#' @examples xxx
gfiUltra <- function(formula, data, gamma = 1, ...){
  y <- f_eval_lhs(formula, data = data)
  X <- model.matrix(formula, data = data)
  n <- length(y)
  p <- ncol(X)
  Xm1 <- X[, -1L, drop = FALSE]
  sis <- SIS(Xm1, y = y, family = "gaussian")
  models <- powerset(sis[["ix"]], colnames(Xm1))
  models_with_FIT <- modelsWithFIT(X = X, y = y, models = models)
  modelsProbs <-
    Rgammas(n = n, p = p, models_with_FIT = models_with_FIT, gamma = gamma)

}
