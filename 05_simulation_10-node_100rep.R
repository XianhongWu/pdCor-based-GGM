# ============================================================
# Script: 05_simulation_10-node_100rep.R
#
# Description:
#   Robust simulation benchmark (p = 10, n = 30, 100 repetitions)
#   comparing two network construction approaches:
#     1) pdCor-based GGM (Partial Distance Correlation)
#     2) Pearson-based GGM (traditional partial correlation)
#
# Workflow:
#   - Generate nonlinear data with known ground-truth edges
#   - Build pdCor and Pearson graphs via significance testing (p < 0.05)
#   - Aggregate selection frequencies across repetitions
#   - Compute PR/ROC curves and AUPRC/AUROC
#   - Plot consensus graphs (≥50% frequency) and heatmaps
#   - Print ranked consensus edges
#
# Dependencies:
#   energy, ppcor, corpcor, qgraph, pheatmap, pROC, PRROC
# ============================================================

### ----------- Load packages ---------------
library(energy)
library(ppcor)
library(corpcor)
library(qgraph)
library(pheatmap)
library(pROC)
library(PRROC)

### ----------- Global settings ---------------
set.seed(123)
n_rep <- 100
n <- 30
p <- 10
vars <- paste0("X", 1:p)
recall_grid <- seq(0, 1, length.out = 100)

# Initialize frequency accumulators
adj_count_pdcor <- matrix(0, p, p, dimnames = list(vars, vars))
adj_count_pearson <- matrix(0, p, p, dimnames = list(vars, vars))

# Initialize storage for PR & ROC curves and AUCs
pr_pdcor_mat <- matrix(0, n_rep, length(recall_grid))
pr_pearson_mat <- matrix(0, n_rep, length(recall_grid))
fpr_grid <- seq(0, 1, length.out = 100)
roc_pdcor_mat <- matrix(NA, n_rep, length(fpr_grid))
roc_pearson_mat <- matrix(NA, n_rep, length(fpr_grid))
auprc_pdcor_all <- numeric(n_rep)
auroc_pdcor_all <- numeric(n_rep)
auprc_pearson_all <- numeric(n_rep)
auroc_pearson_all <- numeric(n_rep)

### ----------- Safety helpers ---------------
safe_cor2pcor <- function(M) {
  if (any(is.na(M)) || any(!is.finite(M))) return(matrix(NA, nrow = nrow(M), ncol = ncol(M)))
  tryCatch(cor2pcor(M), error = function(e) matrix(NA, nrow = nrow(M), ncol = ncol(M)))
}

safe_interp <- function(x, y, grid) {
  if (length(unique(x)) < 2 || any(is.na(x)) || any(is.na(y))) return(rep(NA, length(grid)))
  approx(x, y, xout = grid, rule = 2)$y
}

get_labels_scores <- function(pred_mat, true_mat) {
  idx <- which(upper.tri(true_mat), arr.ind = TRUE)
  y_true <- true_mat[idx]
  y_score <- pred_mat[idx]
  list(y_true = y_true, y_score = y_score)
}

### ----------- Robust simulation loop ---------------
for (rep in 1:n_rep) {
  # Simulate nonlinear data
  noise_sd <- 0.1
  eps <- function(n) rnorm(n, 0, noise_sd)
  X1 <- runif(n, -2, 2)
  X2 <- runif(n, -2, 2)
  X3 <- 4 * X1^2 - 2 + eps(n)
  X4 <- sin(pi * X1) + eps(n)
  X5 <- (2 * X2 + 0.3)^2 - 3 + eps(n)
  X6 <- 0.5 * (X3 - 4)^2 - 1 + eps(n)
  X7 <- 1.2 * log(abs(X3) + 1) + 0.5 * sqrt(abs(X5)) - 3 + eps(n)
  X8 <- 0.3 * X4^3 - 0.4 * (X5 / 4)^2 + 0.8 * X5 + eps(n)
  X9 <- exp(0.2 * (X5 - 3)) - 2 + eps(n)
  X10 <- 0.2 * (X8 / 5)^4 + 1.5 * X2^3 - 4 + eps(n)
  X <- data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)

  # Ground-truth undirected edges
  edges <- data.frame(
    from = c("X1", "X1", "X2", "X3", "X3", "X5", "X5", "X5", "X4", "X8", "X2"),
    to   = c("X3", "X4", "X5", "X6", "X7", "X7", "X8", "X9", "X8", "X10", "X10")
  )
  adj_true <- matrix(0, p, p, dimnames = list(vars, vars))
  for (i in 1:nrow(edges)) {
    adj_true[edges$from[i], edges$to[i]] <- 1
    adj_true[edges$to[i], edges$from[i]] <- 1
  }

  # ----- pdCor path -----
  dCovMat <- matrix(0, p, p)
  for (i in 1:p) for (j in 1:p) dCovMat[i, j] <- dcov(X[[i]], X[[j]])
  dVar <- sqrt(diag(dCovMat))
  dCorMat <- dCovMat / (dVar %*% t(dVar))
  diag(dCorMat) <- 1
  pdCorMat <- safe_cor2pcor(dCorMat)

  pdcor_pval <- matrix(NA, p, p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      z_idx <- setdiff(1:p, c(i, j))
      res <- tryCatch(pdcor.test(X[[i]], X[[j]], as.matrix(X[, z_idx]), R = 500),
                      error = function(e) list(p.value = NA))
      pdcor_pval[i, j] <- pdcor_pval[j, i] <- res$p.value
    }
  }
  adj_tmp_pdcor <- matrix(0, p, p)
  adj_tmp_pdcor[pdcor_pval < 0.05] <- 1
  diag(adj_tmp_pdcor) <- 0
  adj_tmp_pdcor <- (adj_tmp_pdcor + t(adj_tmp_pdcor)) / 2
  adj_count_pdcor <- adj_count_pdcor + adj_tmp_pdcor

  # ----- Pearson path -----
  corMat <- cor(X)
  pcorMat <- safe_cor2pcor(corMat)

  pearson_pval <- matrix(NA, p, p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      z_idx <- setdiff(1:p, c(i, j))
      res <- tryCatch(pcor.test(X[[i]], X[[j]], as.matrix(X[, z_idx])),
                      error = function(e) list(p.value = NA))
      pearson_pval[i, j] <- pearson_pval[j, i] <- res$p.value
    }
  }
  adj_tmp_pearson <- matrix(0, p, p)
  adj_tmp_pearson[pearson_pval < 0.05] <- 1
  diag(adj_tmp_pearson) <- 0
  adj_tmp_pearson <- (adj_tmp_pearson + t(adj_tmp_pearson)) / 2
  adj_count_pearson <- adj_count_pearson + adj_tmp_pearson

  # ----- AUC & PR -----
  pdcor_eval <- get_labels_scores(abs(pdCorMat), adj_true)
  pearson_eval <- get_labels_scores(abs(pcorMat), adj_true)

  if (any(is.na(pdcor_eval$y_score)) || all(is.na(pdcor_eval$y_score))) next
  if (any(is.na(pearson_eval$y_score)) || all(is.na(pearson_eval$y_score))) next

  roc_pdcor <- tryCatch(roc(pdcor_eval$y_true, pdcor_eval$y_score, quiet = TRUE),
                        error = function(e) NULL)
  roc_pearson <- tryCatch(roc(pearson_eval$y_true, pearson_eval$y_score, quiet = TRUE),
                          error = function(e) NULL)

  pr_pdcor <- tryCatch(pr.curve(scores.class0 = pdcor_eval$y_score[pdcor_eval$y_true == 1],
                                scores.class1 = pdcor_eval$y_score[pdcor_eval$y_true == 0],
                                curve = TRUE), error = function(e) NULL)

  pr_pearson <- tryCatch(pr.curve(scores.class0 = pearson_eval$y_score[pearson_eval$y_true == 1],
                                  scores.class1 = pearson_eval$y_score[pearson_eval$y_true == 0],
                                  curve = TRUE), error = function(e) NULL)

  if (!is.null(pr_pdcor)) {
    pr_pdcor_mat[rep, ] <- safe_interp(pr_pdcor$curve[, 1], pr_pdcor$curve[, 2], recall_grid)
    auprc_pdcor_all[rep] <- pr_pdcor$auc.integral
  }
  if (!is.null(pr_pearson)) {
    pr_pearson_mat[rep, ] <- safe_interp(pr_pearson$curve[, 1], pr_pearson$curve[, 2], recall_grid)
    auprc_pearson_all[rep] <- pr_pearson$auc.integral
  }
  if (!is.null(roc_pdcor))  auroc_pdcor_all[rep]  <- auc(roc_pdcor)
  if (!is.null(roc_pearson)) auroc_pearson_all[rep] <- auc(roc_pearson)

  if (!is.null(roc_pdcor)) {
    roc_pdcor_mat[rep, ] <- safe_interp(1 - roc_pdcor$specificities,
                                        roc_pdcor$sensitivities,
                                        fpr_grid)
  }
  if (!is.null(roc_pearson)) {
    roc_pearson_mat[rep, ] <- safe_interp(1 - roc_pearson$specificities,
                                          roc_pearson$sensitivities,
                                          fpr_grid)
  }
}

### ----------- Frequency and consensus graphs ---------------
adj_freq_pdcor <- adj_count_pdcor / n_rep
adj_freq_pearson <- adj_count_pearson / n_rep
adj_consensus_pdcor <- (adj_freq_pdcor >= 0.5) * 1
adj_consensus_pearson <- (adj_freq_pearson >= 0.5) * 1

### ----------- Visualizations ---------------
qgraph(adj_consensus_pdcor, layout = "circle", edge.color = "red",
       labels = vars, title = "Average Graph (pdCor ≥ 50%)")
qgraph(adj_consensus_pearson, layout = "circle", edge.color = "green",
       labels = vars, title = "Average Graph (Pearson ≥ 50%)")

blues <- colorRampPalette(c("#eaf2ff", "#a6c8ff", "#0066ff"))(100)
breaks01 <- seq(0, 1, length.out = 101)

pheatmap(adj_freq_pdcor, cluster_rows = FALSE, cluster_cols = FALSE,
         color = blues, breaks = breaks01, display_numbers = TRUE,
         main = sprintf("Edge Selection Frequency Heatmap (pdCor, n=%d)", n))

pheatmap(adj_freq_pearson, cluster_rows = FALSE, cluster_cols = FALSE,
         color = blues, breaks = breaks01, display_numbers = TRUE,
         main = sprintf("Edge Selection Frequency Heatmap (Pearson, n=%d)", n))

### ----------- Average PR curves ---------------
mean_pdcor_pr <- colMeans(pr_pdcor_mat)
mean_pearson_pr <- colMeans(pr_pearson_mat)

plot(recall_grid, mean_pdcor_pr, type = "l", col = "blue", lwd = 2, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Recall", ylab = "Precision", main = "Average PR Curve (100 runs)")
lines(recall_grid, mean_pearson_pr, col = "red", lwd = 2)
abline(h = mean(pdcor_eval$y_true), col = "gray", lty = 2)
legend("topright", legend = c("pdCor", "Pearson", "Baseline"),
       col = c("blue", "red", "gray"), lwd = 2, lty = c(1, 1, 2))

### ----------- Average ROC curves ---------------
mean_pdcor_roc <- colMeans(roc_pdcor_mat, na.rm = TRUE)
mean_pearson_roc <- colMeans(roc_pearson_mat, na.rm = TRUE)

plot(fpr_grid, mean_pdcor_roc, type = "l", col = "blue", lwd = 2,
     xlim = c(0,1), ylim = c(0,1),
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = "Average ROC Curve (100 runs)")
lines(fpr_grid, mean_pearson_roc, col = "red", lwd = 2)
abline(0, 1, col = "gray", lty = 2)
legend("bottomright", legend = c("pdCor", "Pearson", "Random"),
       col = c("blue", "red", "gray"), lwd = 2, lty = c(1,1,2))

### ----------- Print averaged AUCs ---------------
cat("==== Averaged over", n_rep, "runs ====\n",
    "pdCor  AUPRC:", round(mean(auprc_pdcor_all), 3), " | AUROC:", round(mean(auroc_pdcor_all), 3), "\n",
    "Pearson AUPRC:", round(mean(auprc_pearson_all), 3), " | AUROC:", round(mean(auroc_pearson_all), 3), "\n")

### ----------- Print consensus edges (≥50%), ranked ---------------
print_consensus_edges <- function(freq_mat, method = c("pdCor", "Pearson"),
                                  vars = colnames(freq_mat), cutoff = 0.5, n_rep_total = n_rep) {
  method <- match.arg(method)
  ut <- which(upper.tri(freq_mat), arr.ind = TRUE)
  df <- data.frame(
    method = method,
    node1  = vars[ut[,1]],
    node2  = vars[ut[,2]],
    count  = round(freq_mat[ut] * n_rep_total),
    freq   = freq_mat[ut]
  )
  df <- df[df$freq >= cutoff, ]
  df <- df[order(-df$freq, -df$count), ]
  df$Rank <- seq_len(nrow(df))
  df <- df[, c("Rank", "method", "node1", "node2", "count", "freq")]
  rownames(df) <- NULL
  df
}

cons_pdcor   <- print_consensus_edges(adj_freq_pdcor,  method = "pdCor", vars = vars, cutoff = 0.5)
cons_pearson <- print_consensus_edges(adj_freq_pearson, method = "Pearson", vars = vars, cutoff = 0.5)

cat("\n==== Consensus edges (pdCor, ≥50%) ====\n");    print(cons_pdcor, row.names = FALSE)
cat("\n==== Consensus edges (Pearson, ≥50%) ====\n"); print(cons_pearson, row.names = FALSE)

# (Optional) print all edges ranked by frequency
# 
# print_all_edges <- function(freq_mat, method = c("pdCor", "Pearson"),
#                             vars = colnames(freq_mat), n_rep_total = n_rep) {
#   method <- match.arg(method)
#   ut <- which(upper.tri(freq_mat), arr.ind = TRUE)
#   df <- data.frame(
#     method = method,
#     node1  = vars[ut[,1]],
#     node2  = vars[ut[,2]],
#     count  = round(freq_mat[ut] * n_rep_total),   
#     freq   = freq_mat[ut]                        
#   )
#   df <- df[order(-df$freq, -df$count), ]          
#   df$Rank <- seq_len(nrow(df))                    
#   df <- df[, c("Rank", "method", "node1", "node2", "count", "freq")]
#   rownames(df) <- NULL
#   return(df)
# }
# 
# all_pdcor   <- print_all_edges(adj_freq_pdcor,  method = "pdCor", vars = vars)
# all_pearson <- print_all_edges(adj_freq_pearson, method = "Pearson", vars = vars)
# 
# cat("\n==== All edges ranked by frequency (pdCor) ====\n")
# print(all_pdcor, row.names = FALSE)
# 
# cat("\n==== All edges ranked by frequency (Pearson) ====\n")
# print(all_pearson, row.names = FALSE)
