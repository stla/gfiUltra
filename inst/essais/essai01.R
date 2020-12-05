library(gfiUltra)

# data ####
set.seed(666L)
n <- 10L
x1 <- gl(2L, n/2, labels = c("A", "B"))
x2 <- 1:n
Xnoise <- matrix(rnorm(n*98), nrow = n, ncol = 98)
colnames(Xnoise) <- paste0("x", 3L:100L)
y <- model.matrix(~ x1 + x2) %*% c(10, 20, 3) + rnorm(n, sd = 2)
dat <- cbind(y, data.frame(x1, x2), as.data.frame(Xnoise))

# fiducial simulations ####
fidsims <- gfiUltra(y ~ ., data = dat, nsims = 200)
summary(fidsims)
