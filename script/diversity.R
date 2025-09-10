#!/usr/bin/env Rscript

# Diversity analysis script
# Input: otu.csv and meta.csv files
# Output: Alpha diversity (Chao1, Shannon) and Adonis2 test results

# Load required libraries
library(vegan)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to calculate alpha diversity
calculate_alpha_diversity <- function(otu_data, meta_data) {
  cat("Calculating alpha diversity metrics...\n")
  
  # Ensure samples are in the same order
  common_samples <- intersect(rownames(otu_data), rownames(meta_data))
  otu_data <- otu_data[common_samples, ]
  meta_data <- meta_data[common_samples, ]
  
  # Calculate Chao1 index
  chao1 <- estimateR(otu_data)[2, ]  # Chao1 is the second row
  
  # Calculate Shannon index
  shannon <- diversity(otu_data, index = "shannon")
  
  # Create results dataframe
  alpha_results <- data.frame(
    Sample = rownames(otu_data),
    Chao1 = chao1,
    Shannon = shannon,
    Group = meta_data$group,
    stringsAsFactors = FALSE
  )
  
  return(alpha_results)
}

# Function to perform Adonis2 test
perform_adonis2_test <- function(otu_data, meta_data) {
  cat("Performing Adonis2 test on beta diversity...\n")
  
  # Ensure samples are in the same order
  common_samples <- intersect(rownames(otu_data), rownames(meta_data))
  otu_data <- otu_data[common_samples, ]
  meta_data <- meta_data[common_samples, ]
  
  # Calculate Bray-Curtis distance
  cat("  Calculating Bray-Curtis distance...\n")
  bc_dist <- vegdist(otu_data, method = "bray")
  
  # Calculate Jaccard distance
  cat("  Calculating Jaccard distance...\n")
  jaccard_dist <- vegdist(otu_data, method = "jaccard")
  
  # Perform Adonis2 test on Bray-Curtis distance
  cat("  Performing Adonis2 test on Bray-Curtis distance...\n")
  adonis_bc <- adonis2(bc_dist ~ group, data = meta_data, permutations = 999)
  
  # Perform Adonis2 test on Jaccard distance
  cat("  Performing Adonis2 test on Jaccard distance...\n")
  adonis_jaccard <- adonis2(jaccard_dist ~ group, data = meta_data, permutations = 999)
  
  # Create results list
  adonis_results <- list(
    bray_curtis = adonis_bc,
    jaccard = adonis_jaccard,
    bc_distance_matrix = bc_dist,
    jaccard_distance_matrix = jaccard_dist
  )
  
  return(adonis_results)
}

# Function to create summary statistics
create_summary_stats <- function(alpha_results) {
  cat("Creating summary statistics...\n")
  
  # Summary by group
  summary_stats <- alpha_results %>%
    group_by(Group) %>%
    summarise(
      n_samples = n(),
      Chao1_mean = mean(Chao1, na.rm = TRUE),
      Chao1_sd = sd(Chao1, na.rm = TRUE),
      Chao1_se = sd(Chao1, na.rm = TRUE) / sqrt(n()),
      Shannon_mean = mean(Shannon, na.rm = TRUE),
      Shannon_sd = sd(Shannon, na.rm = TRUE),
      Shannon_se = sd(Shannon, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  return(summary_stats)
}

# Function to perform statistical tests
perform_statistical_tests <- function(alpha_results) {
  cat("Performing statistical tests...\n")
  
  # Test for normality
  chao1_normal <- shapiro.test(alpha_results$Chao1)
  shannon_normal <- shapiro.test(alpha_results$Shannon)
  
  # Perform appropriate statistical test
  if (chao1_normal$p.value > 0.05) {
    chao1_test <- t.test(Chao1 ~ Group, data = alpha_results)
    chao1_method <- "t-test"
  } else {
    chao1_test <- wilcox.test(Chao1 ~ Group, data = alpha_results)
    chao1_method <- "Wilcoxon rank-sum test"
  }
  
  if (shannon_normal$p.value > 0.05) {
    shannon_test <- t.test(Shannon ~ Group, data = alpha_results)
    shannon_method <- "t-test"
  } else {
    shannon_test <- wilcox.test(Shannon ~ Group, data = alpha_results)
    shannon_method <- "Wilcoxon rank-sum test"
  }
  
  # Create results list
  test_results <- list(
    chao1 = list(
      test = chao1_test,
      method = chao1_method,
      normality = chao1_normal
    ),
    shannon = list(
      test = shannon_test,
      method = shannon_method,
      normality = shannon_normal
    )
  )
  
  return(test_results)
}

# Function to create plots
create_plots <- function(alpha_results, output_dir = ".") {
  cat("Creating plots...\n")
  
  # Chao1 boxplot
  p1 <- ggplot(alpha_results, aes(x = Group, y = Chao1, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(title = "Chao1 Index by Group",
         x = "Group",
         y = "Chao1 Index") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Shannon boxplot
  p2 <- ggplot(alpha_results, aes(x = Group, y = Shannon, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(title = "Shannon Index by Group",
         x = "Group",
         y = "Shannon Index") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Save plots
  ggsave(file.path(output_dir, "chao1_boxplot.png"), p1, width = 6, height = 4, dpi = 300)
  ggsave(file.path(output_dir, "shannon_boxplot.png"), p2, width = 6, height = 4, dpi = 300)
  
  cat("  Plots saved to:", output_dir, "\n")
  
  return(list(chao1_plot = p1, shannon_plot = p2))
}

# Main function
main <- function(otu_file, meta_file, output_dir = ".") {
  cat("Starting diversity analysis...\n")
  cat("OTU file:", otu_file, "\n")
  cat("Metadata file:", meta_file, "\n")
  cat("Output directory:", output_dir, "\n\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load data
  cat("Loading data...\n")
  otu_data <- read.csv(otu_file, row.names = 1, check.names = FALSE)
  meta_data <- read.csv(meta_file, row.names = 1, check.names = FALSE)
  
  # Check if group column exists
  if (!"group" %in% colnames(meta_data)) {
    stop("Error: 'group' column not found in metadata file")
  }
  
  # Check if group has normal and sick values
  unique_groups <- unique(meta_data$group)
  if (!all(c("normal", "sick") %in% unique_groups)) {
    cat("Warning: Expected 'normal' and 'sick' groups, but found:", unique_groups, "\n")
  }
  
  cat("  OTU data dimensions:", dim(otu_data), "\n")
  cat("  Metadata dimensions:", dim(meta_data), "\n")
  cat("  Groups found:", unique_groups, "\n\n")
  
  # Calculate alpha diversity
  alpha_results <- calculate_alpha_diversity(otu_data, meta_data)
  
  # Perform Adonis2 test
  adonis_results <- perform_adonis2_test(otu_data, meta_data)
  
  # Create summary statistics
  summary_stats <- create_summary_stats(alpha_results)
  
  # Perform statistical tests
  test_results <- perform_statistical_tests(alpha_results)
  
  # Create plots
  plots <- create_plots(alpha_results, output_dir)
  
  # Save results
  cat("Saving results...\n")
  
  # Save alpha diversity results
  write.csv(alpha_results, file.path(output_dir, "alpha_diversity_results.csv"), row.names = FALSE)
  
  # Save summary statistics
  write.csv(summary_stats, file.path(output_dir, "alpha_diversity_summary.csv"), row.names = FALSE)
  
  # Save Adonis2 results
  sink(file.path(output_dir, "adonis2_results.txt"))
  cat("=== ADONIS2 TEST RESULTS ===\n\n")
  cat("Bray-Curtis Distance:\n")
  print(adonis_results$bray_curtis)
  cat("\nJaccard Distance:\n")
  print(adonis_results$jaccard)
  sink()
  
  # Save statistical test results
  sink(file.path(output_dir, "statistical_tests.txt"))
  cat("=== STATISTICAL TEST RESULTS ===\n\n")
  cat("Chao1 Index:\n")
  cat("Method:", test_results$chao1$method, "\n")
  cat("Normality test p-value:", test_results$chao1$normality$p.value, "\n")
  print(test_results$chao1$test)
  cat("\nShannon Index:\n")
  cat("Method:", test_results$shannon$method, "\n")
  cat("Normality test p-value:", test_results$shannon$normality$p.value, "\n")
  print(test_results$shannon$test)
  sink()
  
  # Print summary to console
  cat("\n=== SUMMARY ===\n")
  cat("Alpha Diversity Summary:\n")
  print(summary_stats)
  
  cat("\nAdonis2 Test Results:\n")
  cat("Bray-Curtis Distance - R²:", adonis_results$bray_curtis$R2[1], 
      ", p-value:", adonis_results$bray_curtis$`Pr(>F)`[1], "\n")
  cat("Jaccard Distance - R²:", adonis_results$jaccard$R2[1], 
      ", p-value:", adonis_results$jaccard$`Pr(>F)`[1], "\n")
  
  cat("\nStatistical Tests:\n")
  cat("Chao1 Index -", test_results$chao1$method, "- p-value:", test_results$chao1$test$p.value, "\n")
  cat("Shannon Index -", test_results$shannon$method, "- p-value:", test_results$shannon$test$p.value, "\n")
  
  cat("\nAnalysis completed! Results saved to:", output_dir, "\n")
  
  return(list(
    alpha_results = alpha_results,
    adonis_results = adonis_results,
    summary_stats = summary_stats,
    test_results = test_results,
    plots = plots
  ))
}

# Command line interface
if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript diversity.R <otu_file> <meta_file> [output_dir]\n")
    cat("Example: Rscript diversity.R otu.csv meta.csv results/\n")
    quit(status = 1)
  }
  
  otu_file <- args[1]
  meta_file <- args[2]
  output_dir <- if (length(args) > 2) args[3] else "."
  
  # Check if files exist
  if (!file.exists(otu_file)) {
    stop("Error: OTU file not found:", otu_file)
  }
  if (!file.exists(meta_file)) {
    stop("Error: Metadata file not found:", meta_file)
  }
  
  # Run analysis
  results <- main(otu_file, meta_file, output_dir)
}

# Example usage (uncomment to use interactively):
# results <- main("otu.csv", "meta.csv", "diversity_results")
