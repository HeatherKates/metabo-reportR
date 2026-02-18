#!/usr/bin/env Rscript

# Set library path
.libPaths(c("R_libs", .libPaths()))

cat("Installing packages to:", .libPaths()[1], "\n\n")

# CRAN packages
cran_packages <- c(
  # Data manipulation
  "dplyr", "tidyr", "data.table", "stringr", "stringi", "plyr",
  
  # File I/O
  "readxl", "xlsx", "yaml", "qs",
  
  # Statistical analysis
  "stats", "car", "broom", "emmeans", "lme4", "lmerTest", "nlme",
  
  # Visualization
  "ggplot2", "plotly", "RColorBrewer", "pheatmap", "heatmaply",
  "gridExtra", "grid",
  
  # Tables and reports
  "DT", "htmltools", "downloadthis", "knitr", "rmarkdown",
  
  # Parallel processing
  "foreach", "parallel", "doParallel",
  
  # Metabolomics-specific
  "omu", "impute"
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
install.packages(cran_packages, lib = .libPaths()[1], repos = "https://cloud.r-project.org/")

# Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = .libPaths()[1])

bioc_packages <- c(
  "MetaboAnalystR",
  "KEGGREST",
  "pathview",
  "org.Hs.eg.db",
  "SummarizedExperiment",
  "S4Vectors",
  "impute"
)

BiocManager::install(bioc_packages, lib = .libPaths()[1], update = FALSE, ask = FALSE)

# GitHub packages
cat("\nInstalling GitHub packages...\n")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", lib = .libPaths()[1])

# Install latest metid from tidymass (replaces local metid_SECIM-main)
remotes::install_github("tidymass/metid", lib = .libPaths()[1], upgrade = "never")

# Quarto (must be installed separately - not an R package)
cat("\n=== NOTE ===\n")
cat("Quarto must be installed separately from: https://quarto.org/docs/get-started/\n")
cat("============\n")

cat("\nPackage installation complete!\n")
cat("Library location:", .libPaths()[1], "\n")
