# 01_introduction.R
# This script generates nonlinear sample data and compares
# Pearson correlation with distance correlation (dCor).
# Used for the Introduction section of the thesis.

# Load required packages
library(energy)  
library(ggplot2)
library(gridExtra)

# Set random seed
set.seed(666)

# Define a function to generate data and plot
plot_non_linear_samples <- function(n) {
  theta <- runif(n, 0, 2 * pi)
  X1 <- cos(theta) + rnorm(n, 0, 0.05)
  X2 <- sin(theta) + rnorm(n, 0, 0.05)

  # Compute correlations
  cor_val <- round(cor(X1, X2), 2)
  dcor_val <- round(dcor(X1, X2), 2)
  
  df <- data.frame(X1 = X1, X2 = X2)
  label_text <- paste0("cor = ", cor_val, "\ndcor = ", dcor_val)
  
  # Plot
  p <- ggplot(df, aes(x = X1, y = X2)) +
    geom_point(shape = 1, size = 2) +
    theme_minimal() +
    labs(x = "X1", y = "X2", title = paste0(n, " samples")) +
    annotate("text", x = min(df$X1), y = max(df$X2), label = label_text, hjust = 0, vjust = 1)
  
  return(p)
}

# Plot with different sample sizes
p1 <- plot_non_linear_samples(1000)
p2 <- plot_non_linear_samples(100)
p3 <- plot_non_linear_samples(20)
p4 <- plot_non_linear_samples(10)

# Arrange plots in a grid
grid.arrange(p1, p2, p3, p4, ncol = 2)
