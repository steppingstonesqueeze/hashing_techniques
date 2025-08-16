
# scripts/run_series_B.R

suppressWarnings(suppressMessages({
  library(optparse)
  library(ggplot2)
  library(Matrix)
  library(matrixStats)
}))

source("R/data_generators.R")
source("R/utils_metrics.R")
source("R/simhash.R")
source("R/pstable_lsh.R")
source("R/plots.R")

option_list <- list(
  make_option(c("--method"), type="character", default="simhash", help="simhash | pstable"),
  make_option(c("--dataset"), type="character", default="iso_gauss", help="iso_gauss | sparse"),
  make_option(c("--n"), type="integer", default=4000, help="number of points"),
  make_option(c("--D"), type="integer", default=256, help="ambient dimension"),
  make_option(c("--pairs"), type="integer", default=80000, help="number of random pairs"),
  make_option(c("--queries"), type="integer", default=200, help="number of query points for ANN eval"),
  make_option(c("--knn_k"), type="integer", default=10, help="k for recall@k"),

  # SimHash params
  make_option(c("--bits"), type="character", default="64,128,256", help="comma-separated bit widths"),

  # p-stable params
  make_option(c("--w"), type="character", default="1.0", help="comma-separated bucket widths"),
  make_option(c("--k_tuple"), type="integer", default=2, help="hashes per table"),
  make_option(c("--L_tables"), type="integer", default=10, help="number of tables"),

  make_option(c("--seed"), type="integer", default=777, help="random seed"),
  make_option(c("--outdir"), type="character", default="results", help="output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$seed)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- data ---
if (opt$dataset == "iso_gauss") {
  X <- gen_gauss_iso(opt$n, opt$D, seed=opt$seed + 1)
} else if (opt$dataset == "sparse") {
  X <- gen_sparse(opt$n, opt$D, density=0.05, seed=opt$seed + 2)
} else {
  stop("unknown dataset")
}

# --- SimHash branch ---
if (opt$method == "simhash") {
  Xu <- normalize_rows(X)  # unit-norm
  pairs <- sample_pairs(nrow(Xu), n_pairs=opt$pairs, seed=opt$seed + 100)
  cos_true <- cosine_similarity_for_pairs(Xu, pairs, chunk=5000)
  d_true <- sqrt(pmax(0, 2 - 2*cos_true))

  bits_list <- as.integer(strsplit(opt$bits, ",")[[1]])
  for (b in bits_list) {
    cat(sprintf("\n[SimHash] bits=%d\n", b))
    B <- simhash_bits(Xu, bits=b, seed=opt$seed + b)
    ham <- hamming_for_pairs(B, pairs, chunk=5000)
    est <- estimates_from_hamming(ham, b)

    # Collision curve (Hamming fraction vs true cosine)
    f_curve <- file.path(opt$outdir, sprintf("simhash_curve_bits%d.png", b))
    plot_collision_curve_simhash(cos_true, est$frac, b, f_curve)

    # Distance estimation diagnostics (unit-norm ℓ2)
    f_scatter <- file.path(opt$outdir, sprintf("simhash_true_vs_est_d2_bits%d.png", b))
    plot_scatter_true_vs_est_d2(d_true, est$d2_hat,
      title=sprintf("SimHash distance estimate (bits=%d)", b),
      outpath=f_scatter)

    relerr <- (est$d2_hat - d_true) / pmax(1e-12, d_true)
    f_hist <- file.path(opt$outdir, sprintf("simhash_relerr_hist_bits%d.png", b))
    plot_hist_error(relerr, title=sprintf("SimHash relative error (bits=%d)", b), outpath=f_hist)

    # Bucket occupancy for grouped keys (bin size effect)
    # Use contiguous blocks of size k_group from the start
    for (k_group in c(8, 16, 32)) {
      if (k_group <= b) {
        keys <- simhash_bucket_keys(B, start_bit=1, k=k_group)
        f_buck <- file.path(opt$outdir, sprintf("simhash_bucket_occupancy_bits%d_kgrp%d.png", b, k_group))
        plot_bucket_occupancy(keys, title=sprintf("SimHash bucket occupancy (bits=%d, key=%d bits)", b, k_group), outpath=f_buck)
      }
    }
  }

  cat("\nSimHash done. See results/ for plots.\n")
  quit(save="no")

} else if (opt$method == "pstable") {
  # --- p-stable branch ---
  pairs <- sample_pairs(nrow(X), n_pairs=opt$pairs, seed=opt$seed + 200)
  d_true <- pairwise_dists_for_pairs(X, pairs, chunk=5000)

  w_list <- as.numeric(strsplit(opt$w, ",")[[1]])
  recall_rows <- list()

  # Build a brute-force KNN baseline for selected queries
  q_idx <- sample.int(nrow(X), opt$queries)
  # Precompute exact distances to all for queries; heavy but ok for n~4000
  exact_nn <- vector("list", length(q_idx))
  for (qi in seq_along(q_idx)) {
    q <- X[q_idx[qi], , drop=FALSE]
    diff <- X - matrix(rep(q, nrow(X)), nrow=nrow(X), byrow=TRUE)
    distv <- sqrt(rowSums(diff*diff))
    ord <- order(distv)
    exact_nn[[qi]] <- ord[2:(opt$knn_k+1)]  # exclude self
  }

  for (w in w_list) {
    cat(sprintf("\n[p-stable] w=%.3f, k=%d, L=%d\n", w, opt$k_tuple, opt$L_tables))
    idx <- build_pstable_index(X, w=w, k_tuple=opt$k_tuple, L_tables=opt$L_tables, seed=opt$seed + as.integer(1000*w))
    # Collision indicator (same bucket in any table) vs true distance
    same_bucket_any <- logical(nrow(pairs))
    for (t in 1:idx$L_tables) {
      A <- idx$A[[t]]; b <- idx$b[[t]]
      Ha <- floor((X[pairs[,1], , drop=FALSE] %*% A + matrix(rep(b, each=nrow(pairs)), nrow=nrow(pairs), byrow=FALSE)) / w)
      Hb <- floor((X[pairs[,2], , drop=FALSE] %*% A + matrix(rep(b, each=nrow(pairs)), nrow=nrow(pairs), byrow=FALSE)) / w)
      same_t <- rowSums(abs(Ha - Hb)) == 0
      same_bucket_any <- same_bucket_any | same_t
    }
    # Plot collision probability vs distance (bin)
    dfc <- data.frame(d=d_true, coll=as.integer(same_bucket_any))
    dfc$bin <- cut(dfc$d, breaks = 50, include.lowest=TRUE)
    agg <- aggregate(coll ~ bin, dfc, mean)
    centers <- sapply(strsplit(as.character(agg$bin), ","), function(s) mean(as.numeric(gsub("\\[|\\]|\\(|\\)", "", s))))
    dplot <- data.frame(center=centers, p=agg$coll)
    p <- ggplot(dplot, aes(x=center, y=p)) + geom_line() + geom_point(size=1) +
      labs(title=sprintf("p-stable: P(same bucket in any table) vs distance (w=%.2f, k=%d, L=%d)", w, opt$k_tuple, opt$L_tables),
           x="true ℓ2 distance", y="collision probability") +
      theme_minimal(base_size=12)
    ggsave(file.path(opt$outdir, sprintf("pstable_collision_vs_d_w%.2f_k%d_L%d.png", w, opt$k_tuple, opt$L_tables)),
           p, width=7, height=5, dpi=140)

    # ANN recall@k
    cand_sizes <- integer(length(q_idx))
    hits <- integer(length(q_idx))
    for (qi in seq_along(q_idx)) {
      id <- q_idx[qi]
      q <- X[id, , drop=FALSE]
      cand <- pstable_candidates(q, idx)
      cand <- setdiff(cand, id)
      cand_sizes[qi] <- length(cand)
      if (length(cand) == 0) {
        hits[qi] <- 0
      } else {
        diff <- X[cand, , drop=FALSE] - matrix(rep(q, length(cand)), nrow=length(cand), byrow=TRUE)
        distv <- sqrt(rowSums(diff*diff))
        ordc <- cand[order(distv)]
        approx_nn <- ordc[1:min(opt$knn_k, length(ordc))]
        hits[qi] <- sum(approx_nn %in% exact_nn[[qi]])
      }
    }
    recall <- mean(hits / opt$knn_k)
    recall_rows[[length(recall_rows)+1]] <- data.frame(method="p-stable", param=sprintf("w=%.2f", w), recall=recall)

    # Candidate size histogram
    plot_candidate_sizes(cand_sizes,
      title=sprintf("p-stable candidate sizes (w=%.2f, k=%d, L=%d)", w, opt$k_tuple, opt$L_tables),
      outpath=file.path(opt$outdir, sprintf("pstable_candidate_sizes_w%.2f_k%d_L%d.png", w, opt$k_tuple, opt$L_tables)))
  }

  dfrec <- do.call(rbind, recall_rows)
  plot_recall_vs_param(dfrec,
    title=sprintf("Recall@%d vs parameter (k=%d, L=%d fixed)", opt$knn_k, opt$k_tuple, opt$L_tables),
    outpath=file.path(opt$outdir, sprintf("pstable_recall_vs_param_k%d_L%d.png", opt$k_tuple, opt$L_tables)))

  cat("\np-stable LSH done. See results/ for plots.\n")
  quit(save="no")

} else {
  stop("unknown method")
}
