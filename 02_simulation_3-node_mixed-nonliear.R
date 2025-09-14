# ============================================================
# Script: 02_simulation_3-node_mixed-nonliear.R
#
# Description:
#   This script simulates a 3-node nonlinear structure and 
#   compares two approaches for conditional dependence testing:
#     1) Partial Distance Correlation (pdCor-based GGM)
#     2) Pearson Partial Correlation (traditional GGM)
#
# Workflow:
#   Step 1. Simulate nonlinear data with 3 variables
#   Step 2. Compute distance correlation matrix and convert to pdCor
#   Step 3. Perform significance testing with pdcor.test()
#   Step 4. Build adjacency matrix for significant pdCor-based network
#   Step 5. Visualize the pdCor network with qgraph
#   Step 6. Fit a traditional Pearson-based GGM for comparison
#   Step 7. Visualize the Pearson network with qgraph
#   Step 8. Compare the two networks using heatmaps
#
# Expected output:
#   - A qgraph plot of the significant pdCor-based network
#   - A qgraph plot of the Pearson-based GGM
#   - Heatmaps highlighting significant edges in both approaches
#
# Dependencies:
#   Packages: energy, corpcor, qgraph, pheatmap, dcov, ppcor
# ============================================================


library(energy)
library(corpcor)
library(qgraph)
library(pheatmap)
library(dcov)
library(ppcor)

### Step 1: Simulate data (nonlinear structure)
set.seed(123)
n <- 300
X1 <- runif(n, -pi, pi)
X2 <- sin(X1) + rnorm(n, 0, 0.1)
X3 <- X2^2 + X1 + rnorm(n, 0, 0.1)
X <- data.frame(X1 = X1, X2 = X2, X3 = X3)
p <- ncol(X)

### Step 2: Construct distance correlation matrix and standardize
dCovMat <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    dCovMat[i, j] <- dcov(X[[i]], X[[j]])
  }
}
dVar <- sqrt(diag(dCovMat))
dCorMat <- dCovMat / (dVar %*% t(dVar))
diag(dCorMat) <- 1

### Step 3: Convert to partial distance correlation matrix
pdCorMat <- cor2pcor(dCorMat)

### Step 4: Build p-value matrix using pdcor.test()
pdCor_pval_matrix <- matrix(NA, p, p)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    z_idx <- setdiff(1:p, c(i, j))
    z <- as.matrix(X[, z_idx])
    res <- pdcor.test(X[[i]], X[[j]], z, R = 500)  # Increase R for more robust results
    pdCor_pval_matrix[i, j] <- pdCor_pval_matrix[j, i] <- res$p.value
  }
}
diag(pdCor_pval_matrix) <- 0

### Step 5: Construct significant network (p < 0.05)
adj_sig <- matrix(0, p, p)
adj_sig[pdCor_pval_matrix < 0.05] <- abs(pdCorMat[pdCor_pval_matrix < 0.05])
diag(adj_sig) <- 0

### Step 6: Visualize significant network
qgraph(adj_sig,
       labels = colnames(X),
       title = "Significant pdCor-based GGM (p < 0.05), n=300",
       fade = FALSE,
       edge.labels = TRUE)

### Step 7: Traditional GGM
pc_res <- ppcor::pcor(as.matrix(X), method = "pearson")
pearson_pcorMat <- pc_res$estimate           # Partial correlation matrix
pearson_pval_matrix <- pc_res$p.value        # Corresponding p-value matrix

# Construct adjacency matrix based on significance (weights = |partial correlation|)
adj_pearson <- matrix(0, p, p)
adj_pearson[pearson_pval_matrix < 0.05] <- abs(pearson_pcorMat[pearson_pval_matrix < 0.05])
diag(adj_pearson) <- 0

# Visualize network
qgraph(adj_pearson,
       labels = colnames(X),
       title = "Pearson-GGM (p < 0.05), n=300",
       fade = FALSE,
       edge.labels = TRUE)

### Step 8: Heatmap comparison
my_palette <- colorRampPalette(c("white", "blue"))(100)
my_breaks <- seq(0, 1, length.out = 101)

# Heatmap of significant pdCor-based GGM
pheatmap(adj_sig,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         main = "Significant pdCor-based GGM (p < 0.05)",
         color = my_palette,
         breaks = my_breaks)

# Heatmap of Pearson-based GGM
pheatmap(adj_pearson,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         main = "Traditional Pearson-GGM (p < 0.05)",
         color = my_palette,
         breaks = my_breaks)
