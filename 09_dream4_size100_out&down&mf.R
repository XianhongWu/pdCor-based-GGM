# ============================================================
# Script: 09_dream4_size100_out&down&mf.R
#
# Description:
#   DREAM4 Size-100 benchmark (5 networks). For each network, this
#   script loads wild-type, knockouts, knockdowns, and multifactorial
#   perturbations, constructs:
#     1) pdCor-based GGM (Partial Distance Correlation)
#     2) Pearson-based GGM (traditional partial correlation),
#   then evaluates both via PR/ROC (AUPRC, AUROC), saves plots and a
#   summary CSV of metrics.
#
# Inputs:
#   - DREAM4 training data folders for Size-100 networks (1..5)
#   - DREAM4 gold standard TSV files (size100_*_goldstandard.tsv)
#   * Paths are set below; adjust to your local machine.
#
# Outputs (per network):
#   - PDFs: True graph, pdCor graph, Pearson graph, PR curve, ROC curve
#   - CSV: dream4_100_out&down&mf_results.csv (summary metrics)
#
# Dependencies:
#   tidyverse, readr, energy, corpcor, qgraph, ppcor, PRROC, pROC
# ============================================================

# ========= Load required libraries ==========
library(tidyverse)
library(readr)
library(energy)
library(corpcor)
library(qgraph)
library(ppcor)
library(PRROC)
library(pROC)

# ========= Function: load gold standard ==========
load_gold_standard <- function(file_path) {
  gold_df <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gold_df) <- c("Source", "Target", "Edge")
  gold_df$Edge <- as.numeric(gold_df$Edge)
  vars <- sort(unique(c(gold_df$Source, gold_df$Target)))
  p <- length(vars)
  adj_mat <- matrix(0, p, p, dimnames = list(vars, vars))
  for (i in 1:nrow(gold_df)) {
    if (gold_df$Edge[i] == 1) {
      from <- gold_df$Source[i]
      to <- gold_df$Target[i]
      adj_mat[from, to] <- 1
      adj_mat[to, from] <- 1
    }
  }
  diag(adj_mat) <- 0
  return(adj_mat)
}

# ========= DREAM4 Size100 batch loop ==========
results_list <- list()

for (net_id in 1:5) {
  
  cat("========== Running Size100 Net", net_id, " (wild+ko+kd+mf) ==========\n")
  
  # ---- 1) Paths ----
  data_path <- paste0("/Users/wuxianhong/Desktop/Thesis/Data-Dream4/Size 100/DREAM4 training data/insilico_size100_", net_id)
  gold_file <- paste0("/Users/wuxianhong/Desktop/Thesis/Data-Dream4/Size 100/DREAM4 gold standards/insilico_size100_", net_id, "_goldstandard.tsv")
  
  # ---- 2) Read and combine data ----
  wildtype <- read_tsv(file.path(data_path, paste0("insilico_size100_", net_id, "_wildtype.tsv")), col_names = TRUE)
  knockouts <- read_tsv(file.path(data_path, paste0("insilico_size100_", net_id, "_knockouts.tsv")), col_names = TRUE)
  knockdowns <- read_tsv(file.path(data_path, paste0("insilico_size100_", net_id, "_knockdowns.tsv")), col_names = TRUE)
  multifactorial <- read_tsv(file.path(data_path, paste0("insilico_size100_", net_id, "_multifactorial.tsv")), col_names = TRUE)
  
  combined_data <- bind_rows(wildtype, knockouts, knockdowns, multifactorial)
  cat("Combined dimensions:", dim(combined_data), "\n")
  
  n <- nrow(combined_data)
  p <- ncol(combined_data)
  
  # ---- 3) Gold standard ----
  adj_true <- load_gold_standard(gold_file)
  adj_true <- adj_true[colnames(combined_data), colnames(combined_data)] # align order
  
  cat("Gold standard edges:", sum(adj_true) / 2, "\n")
  
  pdf(paste0("DREAM4_Size100_Net", net_id, "_out&down&mf_TrueGraph.pdf"), width = 7, height = 7)
  qgraph(adj_true,
         layout = "circle",
         labels = FALSE,
         title = paste0("DREAM4 Size100 Net", net_id, " True Graph(out&down&mf)"),
         fade = FALSE,
         edge.width = 3,
         edge.color = "black",
         vsize = 4,
         color = "white",
         border.color = "black",
         label.color = "black",
         directed = FALSE)
  dev.off()
  
  # ---- 4) pdCor-GGM ----
  cat("Step 1: pdCor calculation...\n")
  dCovMat <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      dCovMat[i, j] <- dcov(combined_data[[i]], combined_data[[j]])
    }
  }
  dVar <- sqrt(diag(dCovMat))
  dCorMat <- dCovMat / (dVar %*% t(dVar))
  diag(dCorMat) <- 1
  pdCorMat <- cor2pcor(dCorMat)
  
  pdCor_pval_matrix <- matrix(NA, p, p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      z_idx <- setdiff(1:p, c(i, j))
      z <- as.matrix(combined_data[, z_idx])
      res <- pdcor.test(combined_data[[i]], combined_data[[j]], z, R = 500)
      pdCor_pval_matrix[i, j] <- pdCor_pval_matrix[j, i] <- res$p.value
    }
  }
  diag(pdCor_pval_matrix) <- 0
  adj_pdcor <- matrix(0, p, p)
  adj_pdcor[pdCor_pval_matrix < 0.05] <- abs(pdCorMat[pdCor_pval_matrix < 0.05])
  diag(adj_pdcor) <- 0
  
  pdf(paste0("DREAM4_Size100_Net", net_id, "_out&down&mf_pdCor.pdf"), width = 7, height = 7)
  qgraph(adj_pdcor,
         layout = "circle",
         labels = FALSE,
         title = paste0("pdCor-GGM Size100 Net", net_id, "(out&down&mf)"),
         fade = FALSE,
         edge.width = TRUE,
         vsize = 4,
         color = "white",
         border.color = "black",
         label.color = "black",
         directed = FALSE)
  dev.off()
  
  # ---- 5) Pearson-GGM ----
  cat("Step 2: Pearson calculation...\n")
  corMat <- cor(combined_data)
  pcorMat <- cor2pcor(corMat)
  
  pearson_pval_matrix <- matrix(NA, p, p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      z_idx <- setdiff(1:p, c(i, j))
      z <- as.matrix(combined_data[, z_idx])
      res <- pcor.test(combined_data[[i]], combined_data[[j]], z)
      pearson_pval_matrix[i, j] <- pearson_pval_matrix[j, i] <- res$p.value
    }
  }
  diag(pearson_pval_matrix) <- 0
  adj_pearson <- matrix(0, p, p)
  adj_pearson[pearson_pval_matrix < 0.05] <- abs(pcorMat[pearson_pval_matrix < 0.05])
  diag(adj_pearson) <- 0
  
  
  pdf(paste0("DREAM4_Size100_Net", net_id, "_out&down&mf_Pearson.pdf"), width = 7, height = 7)
  qgraph(adj_pearson,
         layout = "circle",
         labels = FALSE,
         title = paste0("Pearson-GGM Size100 Net", net_id, "(out&down&mf)"),
         fade = FALSE,
         edge.width = TRUE,
         vsize = 4,
         color = "white",
         border.color = "black",
         label.color = "black",
         directed = FALSE)
  dev.off()
  
  # ---- 6) PR & ROC ----
  cat("Step 3: PR & ROC calculation...\n")
  idx <- upper.tri(adj_true)
  y_true <- adj_true[idx]
  y_pdcor <- adj_pdcor[idx]
  y_pearson <- adj_pearson[idx]
  
  roc_pdcor <- roc(response = y_true, predictor = y_pdcor)
  roc_pearson <- roc(response = y_true, predictor = y_pearson)
  pr_pdcor <- pr.curve(scores.class0 = y_pdcor[y_true == 1], scores.class1 = y_pdcor[y_true == 0], curve = TRUE)
  pr_pearson <- pr.curve(scores.class0 = y_pearson[y_true == 1], scores.class1 = y_pearson[y_true == 0], curve = TRUE)
  
  pdf(paste0("DREAM4_Size100_Net", net_id, "_out&down&mf_PR.pdf"), width = 6.5, height = 5)
  plot(pr_pdcor$curve[,1], pr_pdcor$curve[,2], type = "l", col = "blue", lwd = 2,
       xlab = "Recall", ylab = "Precision", main = paste0("PR Curve Size100 Net", net_id, "(out&down&mf)"))
  lines(pr_pearson$curve[,1], pr_pearson$curve[,2], col = "red", lwd = 2)
  legend("bottomleft", legend = c("pdCor", "Pearson"), col = c("blue", "red"), lty = 1)
  dev.off()
  
  pdf(paste0("DREAM4_Size100_Net", net_id, "_out&down&mf_ROC.pdf"), width = 6, height = 5)
  plot(roc_pdcor, col = "blue", main = paste0("ROC Curve Size100 Net", net_id, "(out&down&mf)"))
  lines(roc_pearson, col = "red")
  legend("bottomright", legend = c("pdCor", "Pearson"), col = c("blue", "red"), lty = 1)
  dev.off()
  
  # ---- 7) Save metrics ----
  results_list[[net_id]] <- tibble(
    Network = paste0("Size100_Net", net_id, "_out&down&mf"),
    Method = c("pdCor", "Pearson"),
    PR_AUC = c(round(pr_pdcor$auc.integral, 4), round(pr_pearson$auc.integral, 4)),
    ROC_AUC = c(round(auc(roc_pdcor), 4), round(auc(roc_pearson), 4))
  )
  
  # Optional: if you want auto-plotting like the 10-node script, add PDF outputs similarly
}

# ========= Aggregate results and save ==========
results_df <- bind_rows(results_list)
write_csv(results_df, "dream4_100_out&down&mf_results.csv")
print(results_df)
