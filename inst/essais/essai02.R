library(gfiUltra)

# data ####
set.seed(666L)
n <- 300L
p <- 1000L
x1 <- gl(2L, n/2, labels = c("A", "B"))
x2 <- 1:n
Xnoise <- matrix(rnorm(n*(p-2L)), nrow = n, ncol = p-2L)
colnames(Xnoise) <- paste0("x", 3L:p)
y <- model.matrix(~ x1 + x2) %*% c(10, 20, 3) + rnorm(n, sd = 0.5)
dat <- cbind(y, data.frame(x1, x2), as.data.frame(Xnoise))

# fiducial simulations ####
gfi <- gfiUltra(y ~ ., data = dat, nsims = 2000)
#summary(gfi$fidSims)

gfiConfInt(gfi)

