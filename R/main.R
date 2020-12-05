#' Generalized fiducial inference for ultrahigh-dimensional regression
#' @description Generates the fiducial simulations of the parameters of
#'   a "large p - small n" regression model and returns the selected models
#'   with their probability.
#'
#' @param formula a formula describing the model
#' @param data dataframe in which to search the variables of the model
#' @param gamma tuning parameter; for expert usage only
#' @param ... named arguments passed to \code{\link[SIS:SIS]{SIS}}, such as
#'   \code{penalty = "lasso"}
#'
#' @return A list with two elements: the fiducial simulations in a matrix and
#'   a vector giving the probabilities of the selected models.
#'
#' @references Randy C. S. Lai, Jan Hannig & Thomas C. M. Lee.
#'   \emph{Generalized Fiducial Inference for Ultrahigh-Dimensional Regression}.
#'   Journal of the American Statistical Association,
#'   Volume 110, 2015 - Issue 510, 760-772.
#'   <doi:10.1080/01621459.2014.931237>
#'
#' @export
#' @importFrom SIS SIS
#' @importFrom lazyeval f_eval_lhs
#' @importFrom stats model.matrix rchisq
#' @importFrom mvtnorm rmvnorm
#'
#' @examples # data ####
#' set.seed(666L)
#' n <- 300L
#' p <- 1000L
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' colnames(X) <- paste0("x", 1L:p)
#' beta <- c(4, 5, 6, 7, 8)
#' y <- X[, 1L:5L] %*% beta + rnorm(n, sd = 0.9)
#' dat <- cbind(y, as.data.frame(X))
#' # fiducial simulations ####
#' gfi <- gfiUltra(y ~ ., data = dat, nsims = 10000L)
#' # fiducial confidence intervals
#' gfiConfInt(gfi)
#' gfiEstimates(gfi)
gfiUltra <- function(formula, data, nsims = 5L, gamma = 1, ...){
  y <- f_eval_lhs(formula, data = data)
  X <- model.matrix(formula, data = data)
  n <- length(y)
  p <- ncol(X)
  Xm1 <- X[, -1L, drop = FALSE]
  sis <- quiet(SIS(Xm1, y = y, family = "gaussian", ...))
  models <- powerset(sis[["ix"]], colnames(Xm1))
  models_with_FIT <- modelsWithFIT(X = X, y = y, models = models)
  modelsProbs <-
    Rgammas(n = n, p = p, models_with_FIT = models_with_FIT, gamma = gamma)
  #
  # simulations
  selected <- colnames(X)[c(1L, 1L + sis[["ix"]])]
  Sims <- matrix(NA_real_, nrow = nsims, ncol = length(selected) + 1L)
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
  list(fidSims = Sims, models = sort(modelsProbs, decreasing = TRUE))
}

#' Fiducial confidence intervals for ultrahigh-dimensional regression
#' @description Fiducial confidence intervals of the selected parameters of a
#'   ultrahigh-dimensional regression.
#'
#' @param gfi an output of \code{\link{gfiUltra}}
#' @param conf confidence level
#'
#' @return The confidence intervals in a matrix.
#' @export
#' @seealso \code{\link{gfiEstimates}}
gfiConfInt <- function(gfi, conf = 0.95){
  Sims <- gfi[["fidSims"]]
  Beta <- Sims[, -ncol(Sims), drop = FALSE]
  NAs <- apply(Beta, 2L, function(x) mean(is.na(x)))
  Beta <- Beta[, NAs < 0.5, drop = FALSE]
  alpha <- 1-conf
  intrvls <- apply(Beta, 2L, function(x){
    quantile(x, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  })
  cbind(
    intrvls,
    sigma = quantile(Sims[, "sigma"], probs = c(alpha/2, 1-alpha/2))
  )
}

#' Fiducial estimates for ultrahigh-dimensional regression
#' @description Fiducial estimates of the selected parameters of a
#'   ultrahigh-dimensional regression.
#'
#' @param gfi an output of \code{\link{gfiUltra}}
#'
#' @return The estimates in a matrix.
#' @export
#' @seealso \code{\link{gfiConfInt}}
gfiEstimates <- function(gfi){
  Sims <- gfi[["fidSims"]]
  Beta <- Sims[, -ncol(Sims), drop = FALSE]
  NAs <- apply(Beta, 2L, function(x) mean(is.na(x)))
  Beta <- Beta[, NAs < 0.5, drop = FALSE]
  estimates <- apply(Beta, 2L, function(x){
    c(median = median(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE))
  })
  cbind(
    estimates,
    sigma = c(median = median(Sims[, "sigma"]), mean = mean(Sims[, "sigma"]))
  )
}
