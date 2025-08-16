
# R/simhash.R

# SimHash (random hyperplanes) for angular/cosine similarity.
# bits: number of hyperplanes (columns of W)
# Returns integer bit matrix (0/1) of size n x bits
simhash_bits <- function(Xu, bits=128, seed=123) {
  set.seed(seed)
  n <- nrow(Xu); D <- ncol(Xu)
  W <- matrix(rnorm(D*bits), D, bits)  # Gaussian hyperplanes
  H <- Xu %*% W
  (H >= 0) * 1L
}

# Hamming distance between two bit rows (integer vectors 0/1)
hamming_for_pairs <- function(B, pairs, chunk=5000) {
  n_pairs <- nrow(pairs)
  out <- numeric(n_pairs)
  start <- 1L
  while (start <= n_pairs) {
    end <- min(n_pairs, start + chunk - 1L)
    idx <- start:end
    a <- B[pairs[idx,1], , drop=FALSE]
    b <- B[pairs[idx,2], , drop=FALSE]
    out[idx] <- rowSums(abs(a - b))
    start <- end + 1L
  }
  out
}

# Estimate angle and cosine distance from Hamming fractions
# frac = H / bits ≈ θ/π  => θ_hat = π * frac
# cos_hat = cos(θ_hat); for unit vectors, d2_hat = sqrt(2 - 2*cos_hat)
estimates_from_hamming <- function(hamming_counts, bits) {
  frac <- hamming_counts / bits
  theta_hat <- pi * frac
  cos_hat <- cos(theta_hat)
  d2_hat <- sqrt(pmax(0, 2 - 2*cos_hat))
  list(frac=frac, theta_hat=theta_hat, cos_hat=cos_hat, d2_hat=d2_hat)
}

# Bucketization by grouping k bits as a key (SimHash "bin size" effect)
# Returns character keys per row for a specific contiguous block of k bits
simhash_bucket_keys <- function(B, start_bit=1, k=16) {
  idx <- start_bit:(start_bit+k-1)
  if (max(idx) > ncol(B)) stop("bit block exceeds bit width")
  sub <- B[, idx, drop=FALSE]
  # Collapse bits to string key
  apply(sub, 1, function(v) paste(v, collapse=""))
}
