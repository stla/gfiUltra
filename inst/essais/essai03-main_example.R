library(gfiUltra)

# data ####
set.seed(666L)
n <- 300L
p <- 1000L
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("x", 1L:p)
beta <- c(4, 5, 6, 7, 8)
y <- X[, 1L:5L] %*% beta + rnorm(n, sd = 0.9)
dat <- cbind(y, as.data.frame(X))
# fiducial simulations ####
gfi <- gfiUltra(y ~ ., data = dat, nsims = 10000L)
# fiducial confidence intervals
gfiConfInt(gfi)

