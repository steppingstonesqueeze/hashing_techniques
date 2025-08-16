
# R/plots.R

suppressWarnings(suppressMessages({
  library(ggplot2)
}))

plot_collision_curve_simhash <- function(true_cos, frac, bits, outpath) {
  # Bin by true cosine, plot mean frac (≈ θ/π) with ribbon (std)
  df <- data.frame(true_cos = true_cos, frac = frac)
  df$bin <- cut(df$true_cos, breaks = seq(-1, 1, by=0.04), include.lowest = TRUE)
  agg <- aggregate(frac ~ bin, df, function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
  mu <- sapply(agg$frac, function(v) v[1]); sd <- sapply(agg$frac, function(v) v[2]); n <- sapply(agg$frac, function(v) v[3])
  centers <- sapply(strsplit(as.character(agg$bin), ","), function(s) mean(as.numeric(gsub("\\[|\\]|\\(|\\)", "", s))))
  out <- data.frame(center = centers, mu = mu, sd = sd, n=n)

  p <- ggplot(out, aes(x=center, y=mu)) +
    geom_line() +
    geom_ribbon(aes(ymin=mu-sd, ymax=mu+sd), alpha=0.2) +
    labs(title=sprintf("SimHash: Hamming fraction vs true cosine (bits=%d)", bits),
         x="true cosine similarity", y="Hamming fraction (≈ angle/π)") +
    theme_minimal(base_size=12)
  ggsave(outpath, p, width=7, height=5, dpi=140)
}

plot_bucket_occupancy <- function(keys, title, outpath) {
  tab <- table(keys)
  df <- data.frame(bucket=names(tab), size=as.vector(tab))
  p <- ggplot(df, aes(x=size)) +
    geom_histogram(bins=60, fill="#2b8cbe", alpha=0.8) +
    labs(title=title, x="bucket size", y="count") +
    theme_minimal(base_size=12)
  ggsave(outpath, p, width=7, height=5, dpi=140)
}

plot_scatter_true_vs_est_d2 <- function(d_true, d_est, title, outpath) {
  k <- length(d_true); take <- if (k > 20000) sample.int(k, 20000) else 1:k
  df <- data.frame(true=d_true[take], est=d_est[take])
  p <- ggplot(df, aes(x=true, y=est)) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    geom_point(alpha=0.2) +
    labs(title=title, x="true ℓ2 distance (unit-norm)", y="estimated from SimHash") +
    theme_minimal(base_size=12)
  ggsave(outpath, p, width=7, height=5, dpi=140)
}

plot_hist_error <- function(err, title, outpath) {
  df <- data.frame(err=err*100.0)
  p <- ggplot(df, aes(x=err)) +
    geom_histogram(bins=60, fill="#7b3294", alpha=0.85) +
    labs(title=title, x="relative error (%)", y="count") +
    theme_minimal(base_size=12)
  ggsave(outpath, p, width=7, height=5, dpi=140)
}

plot_recall_vs_param <- function(df, title, outpath) {
  p <- ggplot(df, aes(x=param, y=recall, group=method, color=method)) +
    geom_line() + geom_point() +
    labs(title=title, x="parameter", y="recall@k") +
    theme_minimal(base_size=12)
  ggsave(outpath, p, width=7, height=5, dpi=140)
}

plot_candidate_sizes <- function(sizes, title, outpath) {
  df <- data.frame(size=sizes)
  p <- ggplot(df, aes(x=size)) +
    geom_histogram(bins=60, fill="#238b45", alpha=0.85) +
    labs(title=title, x="candidate set size", y="count") +
    theme_minimal(base_size=12)
  ggsave(outpath, p, width=7, height=5, dpi=140)
}
