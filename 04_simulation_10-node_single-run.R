# ============================================================
# Script: 04_simulation_10-node_single-run.R
#
# Description:
#   Single-run simulation (n = 30, p = 10) with nonlinear relations.
#   Builds two graphs via significance testing:
#     1) pdCor-based GGM (Partial Distance Correlation)
#     2) Pearson-based GGM (traditional partial correlation)
#   Outputs:
#     - qgraph visualizations for both methods
#     - edge lists of significant edges (sorted by |weight|)
#     - heatmaps of significant adjacency matrices
#
# Dependencies:
#   energy, corpcor, qgraph, pheatmap, ppcor
# ============================================================

### ----------- Load packages ---------------
library(energy)
library(corpcor)
library(qgraph)
library(pheatmap)
library(ppcor)

### ----------- Step 1: Simulate 10D nonlinear data ---------------
set.seed(123)
n <- 30
noise_sd <- 0.1
eps <- function(n) rnorm(n, 0, noise_sd)

X1  <- runif(n, -2, 2) 
X2  <- runif(n, -2, 2)
X3  <- 4 * X1^2 - 2 + eps(n)                         # quadratic
X4  <- sin(pi * X1) + eps(n)                         # sine
X5  <- (2 * X2 + 0.3)^2 - 3 + eps(n)                 # quadratic + shift
X6  <- 0.5 * (X3 - 4)^2 - 1 + eps(n)                 # chain quadratic
X7  <- 1.2 * log(abs(X3) + 1) + 0.5 * sqrt(abs(X5)) - 3 + eps(n)  # collider: log + sqrt
X8  <- 0.3 * X4^3 - 0.4 * (X5 / 4)^2 + 0.8 * X5 + eps(n)          # collider: cubic + quadratic + linear
X9  <- exp(0.2 * (X5 - 3)) - 2 + eps(n)              # exponential
X10 <- 0.2 * (X8 / 5)^4 + 1.5 * X2^3 - 4 + eps(n)    # quartic + cubic

X <- data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
p <- ncol(X)

### ----------- Step 2: pdCor-GGM construction ---------------
# Build distance correlation matrix
dCovMat <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    dCovMat[i, j] <- dcov(X[[i]], X[[j]])
  }
}
dVar <- sqrt(diag(dCovMat))
dCorMat <- dCovMat / (dVar %*% t(dVar))
diag(dCorMat) <- 1

# Convert to partial distance correlation
pdCorMat <- cor2pcor(dCorMat)

# p-value matrix via pdcor.test()
pdCor_pval_matrix <- matrix(NA, p, p)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    z_idx <- setdiff(1:p, c(i, j))
    z <- as.matrix(X[, z_idx])
    res <- pdcor.test(X[[i]], X[[j]], z, R = 1000)
    pdCor_pval_matrix[i, j] <- pdCor_pval_matrix[j, i] <- res$p.value
  }
}
diag(pdCor_pval_matrix) <- 0

# Symmetrize to avoid double-directed edges in qgraph
pdCorMat <- (pdCorMat + t(pdCorMat)) / 2

# Build adjacency matrix of significant edges
adj_pdcor <- matrix(0, p, p)
adj_pdcor[pdCor_pval_matrix < 0.05] <- abs(pdCorMat[pdCor_pval_matrix < 0.05])
diag(adj_pdcor) <- 0

# Visualize pdCor network
qgraph(adj_pdcor,
       labels = colnames(X),
       title = "Significant pdCor-based GGM (p < 0.05), n=30, p=10",
       fade = FALSE,
       edge.labels = TRUE,
       directed = FALSE)

# ---------- Export pdCor significant edge list ----------
# Use upper triangle only to avoid duplicates
idx <- which(pdCor_pval_matrix < 0.05 & upper.tri(pdCor_pval_matrix), arr.ind = TRUE)

edge_list_pdcor <- data.frame(
  node1   = colnames(X)[idx[, 1]],
  node2   = colnames(X)[idx[, 2]],
  weight  = abs(pdCorMat)[idx],
  p_value = pdCor_pval_matrix[idx]
)

# Sort by absolute weight (descending)
edge_list_pdcor <- edge_list_pdcor[order(edge_list_pdcor$weight, decreasing = TRUE), ]
rownames(edge_list_pdcor) <- NULL
print(edge_list_pdcor)

### ----------- Step 3: Pearson-GGM construction ---------------
corMat <- cor(X)
pcorMat <- cor2pcor(corMat)

# p-value matrix via pcor.test()
pearson_pval_matrix <- matrix(NA, p, p)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    z_idx <- setdiff(1:p, c(i, j))
    z <- as.matrix(X[, z_idx])
    res <- pcor.test(X[[i]], X[[j]], z)
    pearson_pval_matrix[i, j] <- pearson_pval_matrix[j, i] <- res$p.value
  }
}
diag(pearson_pval_matrix) <- 0

# Symmetrize, then build adjacency of significant edges
pcorMat <- (pcorMat + t(pcorMat)) / 2
adj_pearson <- matrix(0, p, p)
adj_pearson[pearson_pval_matrix < 0.05] <- abs(pcorMat[pearson_pval_matrix < 0.05])
diag(adj_pearson) <- 0

# Visualize Pearson network
qgraph(adj_pearson,
       labels = colnames(X),
       title = "Significant Pearson-based GGM (p < 0.05), n=30, p=10",
       fade = FALSE,
       edge.labels = TRUE,
       directed = FALSE)

# ---------- Export Pearson significant edge list ----------
idx <- which(pearson_pval_matrix < 0.05 & upper.tri(pearson_pval_matrix), arr.ind = TRUE)

edge_list_pearson <- data.frame(
  node1   = colnames(X)[idx[, 1]],
  node2   = colnames(X)[idx[, 2]],
  weight  = abs(pcorMat)[idx],
  p_value = pearson_pval_matrix[idx]
)

edge_list_pearson <- edge_list_pearson[order(edge_list_pearson$weight, decreasing = TRUE), ]
rownames(edge_list_pearson) <- NULL
print(edge_list_pearson)

### ----------- Step 4: Heatmap comparison -------------------
# Common palette and breaks
my_palette <- colorRampPalette(c("white", "blue"))(100)
my_breaks <- seq(0, 1, length.out = 101)

# Heatmap for pdCor significant adjacency
pheatmap(adj_pdcor,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         main = "Significant pdCor-based GGM (p < 0.05), n=30",
         color = my_palette,
         breaks = my_breaks,
         fontsize_row = 8,
         fontsize_col = 8)

# Heatmap for Pearson significant adjacency
pheatmap(adj_pearson,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         main = "Significant Pearson-based GGM (p < 0.05), n=30",
         color = my_palette,
         breaks = my_breaks,
         fontsize_row = 8,
         fontsize_col = 8)
