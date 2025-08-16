
# R/utils_metrics.R

suppressWarnings(suppressMessages({
  library(matrixStats)
}))

sample_pairs <- function(n, n_pairs, seed=1234) {
  set.seed(seed)
  i <- sample.int(n, n_pairs, replace=TRUE)
  j <- sample.int(n, n_pairs, replace=TRUE)
  mask <- which(i != j)
  cbind(i[mask], j[mask])
}

pairwise_dists_for_pairs <- function(X, pairs, chunk=5000) {
  n_pairs <- nrow(pairs)
  out <- numeric(n_pairs)
  start <- 1L
  while (start <= n_pairs) {
    end <- min(n_pairs, start + chunk - 1L)
    idx <- start:end
    diff <- X[pairs[idx,1], , drop=FALSE] - X[pairs[idx,2], , drop=FALSE]
    out[idx] <- sqrt(rowSums(diff * diff))
    start <- end + 1L
  }
  out
}

cosine_similarity_for_pairs <- function(Xu, pairs, chunk=5000) {
  # Xu must be unit-normalized
  n_pairs <- nrow(pairs)
  out <- numeric(n_pairs)
  start <- 1L
  while (start <= n_pairs) {
    end <- min(n_pairs, start + chunk - 1L)
    idx <- start:end
    a <- Xu[pairs[idx,1], , drop=FALSE]
    b <- Xu[pairs[idx,2], , drop=FALSE]
    out[idx] <- rowSums(a*b) # dot since unit-norm
    start <- end + 1L
  }
  out
}
