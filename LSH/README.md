
# Series B — LSH in R (SimHash for cosine / angular, p-stable LSH for ℓ2)

This repo demonstrates similarity hashing **(not dimensionality reduction)** and produces:
- Collision probability vs true similarity curves
- Before/after *estimated* distance error (from sketches) as scatter and histograms
- Bucket occupancy histograms and candidate set diagnostics
- ANN recall@k vs parameters (bits/tables for SimHash; (w,k,L) for p-stable)

## Quick start

SimHash (cosine / angular):
```bash
Rscript scripts/run_series_B.R --method simhash --dataset iso_gauss --n 4000 --D 256 \
  --bits 32,64,128,256 --pairs 80000 --queries 200 --knn_k 10
```

p-stable LSH (ℓ2):
```bash
Rscript scripts/run_series_B.R --method pstable --dataset iso_gauss --n 4000 --D 256 \
  --w 0.5,1.0,2.0 --k_tuple 2 --L_tables 10 --pairs 80000 --queries 200 --knn_k 10
```

Outputs (PNGs) are written to `results/`.

### Dependencies
- base R (≥ 4.0)
- packages: `ggplot2`, `Matrix`, `matrixStats`, `optparse`
  ```r
  install.packages(c("ggplot2", "Matrix", "matrixStats", "optparse"))
  ```

### Notes
- For **SimHash**, we **unit-normalize rows** because angular/cosine similarity is the target.
- For **p-stable**, we do not normalize; it targets Euclidean.
- Pairwise curves are computed on a sampled set of pairs for scalability.
- ANN recall@k uses a subset of queries and measures candidate sizes and recall.
