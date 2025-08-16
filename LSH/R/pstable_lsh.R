
# R/pstable_lsh.R

# p-stable LSH for ℓ2 (Datar et al.)
# Hash: h(x) = floor((a · x + b)/w), where a ~ N(0, I), b ~ U(0,w)
# We use k such hashes per table; L tables.
# Index stores map from tuple key to vector indices.

build_pstable_index <- function(X, w=1.0, k_tuple=2, L_tables=10, seed=321) {
  set.seed(seed)
  n <- nrow(X); D <- ncol(X)
  A_list <- vector("list", L_tables)
  b_list <- vector("list", L_tables)
  tables <- vector("list", L_tables)

  for (t in 1:L_tables) {
    A <- matrix(rnorm(D * k_tuple), D, k_tuple)
    b <- runif(k_tuple, min=0, max=w)
    # compute k hashes
    H <- floor((X %*% A + matrix(rep(b, each=n), n, k_tuple, byrow=FALSE)) / w)
    # build keys as "h1|h2|...|hk"
    keys <- apply(H, 1, function(v) paste(v, collapse="|"))
    # bucket
    tab <- split(seq_len(n), keys)
    A_list[[t]] <- A
    b_list[[t]] <- b
    tables[[t]] <- tab
  }
  list(A=A_list, b=b_list, tables=tables, w=w, k_tuple=k_tuple, L_tables=L_tables)
}

# Retrieve candidate set for a query q across all tables
pstable_candidates <- function(q, index) {
  cand <- integer(0)
  for (t in 1:index$L_tables) {
    A <- index$A[[t]]
    b <- index$b[[t]]
    Hq <- floor(((as.numeric(q %*% A)) + b) / index$w)
    key <- paste(Hq, collapse="|")
    bucket <- index$tables[[t]][[key]]
    if (!is.null(bucket)) cand <- c(cand, bucket)
  }
  unique(cand)
}
