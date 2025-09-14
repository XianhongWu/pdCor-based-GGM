# ============================================================
# Script: 03_simulation_3-node_circular.R
#
# Description:
#   This script simulates a 3-node circular structure and 
#   compares two approaches for network construction:
#     1) Partial Distance Correlation (pdCor-based GGM)
#     2) Pearson Partial Correlation (traditional GGM)
#
# Workflow:
#   Step 1. Simulate circular structure data (X1, X2, X3)
#   Step 2. Compute distance correlation matrix (dCor)
#   Step 3. Convert to partial distance correlation matrix (pdCor)
#   Step 4. Perform significance testing with pdcor.test()
#   Step 5. Construct adjacency matrix for significant pdCor-based network
#   Step 6. Visualize pdCor-based network using qgraph
#   Step 7. Construct and visualize traditional Pearson-based GGM
#
# Expected output:
#   - A qgraph plot of the significant pdCor-based network
#   - A qgraph plot of the Pearson-based GGM
#
# Dependencies:
#   Packages: energy, corpcor, qgraph, dcov, ppcor
# ============================================================

### ----------- Simulate data (circular structure) ---------------
set.seed(42)
n <- 300
theta <- runif(n, 0, 2 * pi)
X1 <- cos(theta)
X2 <- sin(theta)
X3 <- X1^2 - X2^2 + rnorm(n, 0, 0.1)
X <- data.frame(X1 = X1, X2 = X2, X3 = X3)
p <- ncol(X)

### ----------- Load required packages ---------------
library(energy)
library(corpcor)
library(qgraph)
library(dcov)

### ----------- Distance Correlation matrix ---------------
dCovMat <- matrix(0, p, p)
for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    dCovMat[i,j] <- dCovMat[j,i] <- dcov(X[[i]], X[[j]])
  }
}
for (i in 1:p) {
  dCovMat[i, i] <- dcov(X[[i]], X[[i]])
}

dVar <- sqrt(diag(dCovMat))
dCorMat <- dCovMat / (dVar %*% t(dVar))
diag(dCorMat) <- 1

### ----------- Partial Distance Correlation matrix ---------------
pdCorMat <- cor2pcor(dCorMat)

### ----------- Significance testing with pdcor.test() ---------------
pdCor_pval_matrix <- matrix(NA, p, p)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    z_idx <- setdiff(1:p, c(i, j))
    z <- as.matrix(X[, z_idx])
    res <- pdcor.test(X[[i]], X[[j]], z, R = 500)  # Increase R for more stable results
    pdCor_pval_matrix[i, j] <- pdCor_pval_matrix[j, i] <- res$p.value
  }
}
diag(pdCor_pval_matrix) <- 0

### ----------- Construct significant network (p < 0.05) ---------------
adj_sig <- matrix(0, p, p)
adj_sig[pdCor_pval_matrix < 0.05] <- abs(pdCorMat[pdCor_pval_matrix < 0.05])
diag(adj_sig) <- 0

### ----------- Visualize significant pdCor-based network ---------------
qgraph(adj_sig,
       labels = colnames(X),
       title = "Circular Structure: Significant pdCor-based GGM (p < 0.05), n=300",
       fade = FALSE,
       edge.labels = TRUE)

### ----------- Traditional GGM (Pearson partial correlation) ---------------
pc_res <- ppcor::pcor(as.matrix(X), method = "pearson")
pearson_pcorMat <- pc_res$estimate       
pearson_pval_matrix <- pc_res$p.value     

# Construct adjacency matrix based on significance (weights = |partial correlation|)
adj_pearson <- matrix(0, p, p)
adj_pearson[pearson_pval_matrix < 0.05] <- abs(pearson_pcorMat[pearson_pval_matrix < 0.05])
diag(adj_pearson) <- 0

### ----------- Visualize Pearson-based GGM ---------------
qgraph(adj_pearson,
       labels = colnames(X),
       title = "Circular Structure: Pearson-GGM (p < 0.05), n=300",
       fade = FALSE,
       edge.labels = TRUE)
