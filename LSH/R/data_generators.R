
# R/data_generators.R

normalize_rows <- function(X) {
  nr <- sqrt(rowSums(X*X))
  nr[nr == 0] <- 1
  X / nr
}

gen_gauss_iso <- function(n, D, seed=42) {
  set.seed(seed)
  matrix(rnorm(n*D), nrow=n, ncol=D)
}

gen_sparse <- function(n, D, density=0.05, seed=44) {
  set.seed(seed)
  X <- matrix(0, n, D)
  nnz <- max(1, floor(density * D))
  for (i in 1:n) {
    idx <- sample.int(D, nnz)
    X[i, idx] <- rnorm(nnz)
  }
  X
}
