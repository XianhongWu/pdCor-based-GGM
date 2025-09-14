# pdCor-based-GGM
R scripts and code to reproduce the experiments in the Master's thesis: "On Distance Graphical Models: A new graphical modeling framework for complex data using partial distance correlations"

The project implements and compares Partial Distance Correlation (pdCor)-based Gaussian Graphical Models (GGMs) with traditional Pearson-based GGMs, through simulations and benchmark datasets (DREAM4 challenge).

⸻
File Organization

	•	01_* → Introduction (demonstration of distance correlation)
	•	02–05_* → Simulation studies (3-node mixed-nonlinear/circular structure, 10-node examples, robustness checks, etc.)
	•	06–09_* → Applications (DREAM4 challenge networks for size10/size100)

⸻
Requirements

R (≥ 4.2) with the following packages:

install.packages(c(
  "energy", "corpcor", "ppcor", "qgraph", "pheatmap",
  "pROC", "PRROC", "tidyverse", "readr", "gridExtra",
  "ggplot2", "dcov"
))


⸻
Usage

1.	Copy or download this repository.
2.	Open R or RStudio, and set the working directory to the project folder.
3.	Run each script independently, depending on your purpose.

