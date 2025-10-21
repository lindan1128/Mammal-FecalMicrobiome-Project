#!/usr/bin/env Rscript

# Differential abundance analysis script for meta-analysis
# Input: 3 OTU files and 3 metadata files
# Output: Differential microbes based on 3 selection criteria

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(reshape2)

# Function to load and process individual dataset
load_dataset <- function(otu_file, meta_file, dataset_name) {
  cat("Loading dataset:", dataset_name, "\n")
  
  # Load OTU data
  otu_data <- read.csv(otu_file, row.names = 1, check.names = FALSE)
  cat("  OTU data dimensions:", dim(otu_data), "\n")
  
  # Load metadata
  meta_data <- read.csv(meta_file, row.names = 1, check.names = FALSE)
  cat("  Metadata dimensions:", dim(meta_data), "\n")
  
  # Check if group column exists
  if (!"group" %in% colnames(meta_data)) {
    stop("Error: 'group' column not found in metadata file:", meta_file)
  }
  
  # Check if group has normal and sick values
  unique_groups <- unique(meta_data$group)
  cat("  Groups found:", unique_groups, "\n")
  
  # Find common samples
  common_samples <- intersect(colnames(otu_data), rownames(meta_data))
  otu_data <- otu_data[, common_samples]
  meta_data <- meta_data[common_samples, ]
  
  cat("  Common samples:", length(common_samples), "\n")
  
  return(list(
    otu_data = otu_data,
    meta_data = meta_data,
    dataset_name = dataset_name
  ))
}

# Function to merge multiple datasets
merge_datasets <- function(datasets) {
  cat("\nMerging datasets...\n")
  
  # Get union of all taxa (row names)
  all_taxa <- c()
  for (dataset in datasets) {
    all_taxa <- union(all_taxa, rownames(dataset$otu_data))
  }
  cat("Total unique taxa across all datasets:", length(all_taxa), "\n")
  
  # Merge OTU data
  merged_otu <- data.frame(row.names = all_taxa)
  merged_meta <- data.frame()
  
  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    dataset_name <- dataset$dataset_name
    
    # Add dataset prefix to sample names to avoid conflicts
    new_sample_names <- paste0(dataset_name, "_", colnames(dataset$otu_data))
    colnames(dataset$otu_data) <- new_sample_names
    rownames(dataset$meta_data) <- new_sample_names
    
    # Add dataset column to metadata
    dataset$meta_data$dataset <- dataset_name
    
    # Merge OTU data (fill missing taxa with 0)
    dataset_otu_aligned <- dataset$otu_data[all_taxa, , drop = FALSE]
    dataset_otu_aligned[is.na(dataset_otu_aligned)] <- 0
    rownames(dataset_otu_aligned) <- all_taxa
    
    merged_otu <- cbind(merged_otu, dataset_otu_aligned)
    merged_meta <- rbind(merged_meta, dataset$meta_data)
  }
  
  cat("Merged OTU data dimensions:", dim(merged_otu), "\n")
  cat("Merged metadata dimensions:", dim(merged_meta), "\n")
  cat("Groups in merged data:", table(merged_meta$group), "\n")
  
  return(list(
    otu_data = merged_otu,
    meta_data = merged_meta
  ))
}

# Function to use relative abundance data (already converted)
use_relative_abundance <- function(otu_data) {
  cat("Using relative abundance data (already converted)...\n")
  
  # Data is already in relative abundance format
  rel_abundance <- otu_data
  
  return(rel_abundance)
}

# Function to apply filtering criteria
apply_filtering_criteria <- function(rel_abundance, meta_data) {
  cat("Applying filtering criteria...\n")
  
  # Criterion 1: Remove taxa with relative abundance < 0.001 across all samples
  cat("  Criterion 1: Filtering taxa with relative abundance < 0.001...\n")
  max_rel_abundance <- apply(rel_abundance, 1, max)
  taxa_to_keep <- max_rel_abundance >= 0.001
  
  filtered_rel_abundance <- rel_abundance[taxa_to_keep, , drop = FALSE]
  cat("  Taxa before filtering:", nrow(rel_abundance), "\n")
  cat("  Taxa after filtering:", nrow(filtered_rel_abundance), "\n")
  
  return(filtered_rel_abundance)
}

# Function to calculate log2 fold change
calculate_log2_fold_change <- function(rel_abundance, meta_data) {
  cat("Calculating log2 fold change...\n")
  
  # Separate cases and controls
  case_samples <- rownames(meta_data)[meta_data$group == "sick"]
  control_samples <- rownames(meta_data)[meta_data$group == "normal"]
  
  cat("  Case samples:", length(case_samples), "\n")
  cat("  Control samples:", length(control_samples), "\n")
  
  # Calculate mean relative abundance for each group
  case_means <- rowMeans(rel_abundance[, case_samples, drop = FALSE], na.rm = TRUE)
  control_means <- rowMeans(rel_abundance[, control_samples, drop = FALSE], na.rm = TRUE)
  
  # Add small pseudocount to avoid log(0)
  pseudocount <- 1e-6
  case_means <- case_means + pseudocount
  control_means <- control_means + pseudocount
  
  # Calculate log2 fold change (case/control)
  log2_fc <- log2(case_means / control_means)
  
  return(log2_fc)
}

# Function to perform Wilcoxon rank sum test
perform_wilcoxon_test <- function(rel_abundance, meta_data) {
  cat("Performing Wilcoxon rank sum test...\n")
  
  # Separate cases and controls
  case_samples <- rownames(meta_data)[meta_data$group == "sick"]
  control_samples <- rownames(meta_data)[meta_data$group == "normal"]
  
  # Initialize results
  p_values <- numeric(nrow(rel_abundance))
  names(p_values) <- rownames(rel_abundance)
  
  # Perform test for each taxon
  for (i in 1:nrow(rel_abundance)) {
    taxon_name <- rownames(rel_abundance)[i]
    case_values <- as.numeric(rel_abundance[i, case_samples])
    control_values <- as.numeric(rel_abundance[i, control_samples])
    
    # Perform Wilcoxon rank sum test
    test_result <- wilcox.test(case_values, control_values, alternative = "two.sided")
    p_values[i] <- test_result$p.value
  }
  
  # Apply Benjamini-Hochberg correction
  cat("  Applying Benjamini-Hochberg correction...\n")
  q_values <- p.adjust(p_values, method = "BH")
  
  return(list(
    p_values = p_values,
    q_values = q_values
  ))
}

# Function to identify differential microbes
identify_differential_microbes <- function(rel_abundance, meta_data) {
  cat("Identifying differential microbes...\n")
  
  # Apply filtering criteria
  filtered_rel_abundance <- apply_filtering_criteria(rel_abundance, meta_data)
  
  # Calculate log2 fold change
  log2_fc <- calculate_log2_fold_change(filtered_rel_abundance, meta_data)
  
  # Criterion 2: Filter by log2 fold change
  cat("  Criterion 2: Filtering by log2 fold change (|log2FC| >= 1)...\n")
  fc_filter <- abs(log2_fc) >= 1
  fc_filtered_rel_abundance <- filtered_rel_abundance[fc_filter, , drop = FALSE]
  fc_filtered_log2_fc <- log2_fc[fc_filter]
  
  cat("  Taxa after log2FC filtering:", nrow(fc_filtered_rel_abundance), "\n")
  
  # Perform statistical test
  test_results <- perform_wilcoxon_test(fc_filtered_rel_abundance, meta_data)
  
  # Criterion 3: Filter by adjusted p-value
  cat("  Criterion 3: Filtering by adjusted p-value (q < 0.05)...\n")
  significant_taxa <- test_results$q_values < 0.05
  significant_rel_abundance <- fc_filtered_rel_abundance[significant_taxa, , drop = FALSE]
  significant_log2_fc <- fc_filtered_log2_fc[significant_taxa]
  significant_p_values <- test_results$p_values[significant_taxa]
  significant_q_values <- test_results$q_values[significant_taxa]
  
  cat("  Final differential microbes:", length(significant_taxa), "\n")
  
  # Create results dataframe
  differential_results <- data.frame(
    taxon = rownames(significant_rel_abundance),
    log2_fold_change = significant_log2_fc,
    p_value = significant_p_values,
    q_value = significant_q_values,
    stringsAsFactors = FALSE
  )
  
  # Add mean relative abundance for each group
  case_samples <- rownames(meta_data)[meta_data$group == "sick"]
  control_samples <- rownames(meta_data)[meta_data$group == "normal"]
  
  case_means <- rowMeans(significant_rel_abundance[, case_samples, drop = FALSE], na.rm = TRUE)
  control_means <- rowMeans(significant_rel_abundance[, control_samples, drop = FALSE], na.rm = TRUE)
  
  differential_results$case_mean_abundance <- case_means
  differential_results$control_mean_abundance <- control_means
  
  # Sort by absolute log2 fold change
  differential_results <- differential_results[order(abs(differential_results$log2_fold_change), decreasing = TRUE), ]
  
  return(differential_results)
}

# Function to create summary statistics
create_summary_stats <- function(differential_results) {
  cat("Creating summary statistics...\n")
  
  # Count up-regulated and down-regulated microbes
  up_regulated <- sum(differential_results$log2_fold_change > 0)
  down_regulated <- sum(differential_results$log2_fold_change < 0)
  
  # Summary statistics
  summary_stats <- list(
    total_differential_microbes = nrow(differential_results),
    up_regulated = up_regulated,
    down_regulated = down_regulated,
    mean_log2_fc = mean(differential_results$log2_fold_change),
    median_log2_fc = median(differential_results$log2_fold_change),
    min_q_value = min(differential_results$q_value),
    max_q_value = max(differential_results$q_value)
  )
  
  return(summary_stats)
}

# Function to create visualizations
create_visualizations <- function(differential_results, output_dir) {
  cat("Creating visualizations...\n")
  
  # Volcano plot
  p1 <- ggplot(differential_results, aes(x = log2_fold_change, y = -log10(q_value))) +
    geom_point(aes(color = ifelse(q_value < 0.05, "Significant", "Not significant")), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    labs(title = "Volcano Plot of Differential Microbes",
         x = "Log2 Fold Change (Case/Control)",
         y = "-Log10(Q-value)",
         color = "Significance") +
    theme_minimal() +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray"))
  
  # Top differential microbes bar plot
  top_microbes <- head(differential_results, 20)
  top_microbes$taxon <- factor(top_microbes$taxon, levels = rev(top_microbes$taxon))
  
  p2 <- ggplot(top_microbes, aes(x = taxon, y = log2_fold_change, 
                                 fill = ifelse(log2_fold_change > 0, "Up-regulated", "Down-regulated"))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Top 20 Differential Microbes",
         x = "Taxon",
         y = "Log2 Fold Change (Case/Control)",
         fill = "Regulation") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +
    scale_fill_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue"))
  
  # Save plots
  ggsave(file.path(output_dir, "volcano_plot.png"), p1, width = 10, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "top_differential_microbes.png"), p2, width = 10, height = 8, dpi = 300)
  
  cat("  Plots saved to:", output_dir, "\n")
  
  return(list(volcano_plot = p1, top_microbes_plot = p2))
}

# Main function
main <- function(otu_files, meta_files, dataset_names, output_dir = ".") {
  cat("Starting differential abundance analysis...\n")
  cat("Number of datasets:", length(otu_files), "\n")
  cat("Output directory:", output_dir, "\n\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load all datasets
  datasets <- list()
  for (i in seq_along(otu_files)) {
    datasets[[i]] <- load_dataset(otu_files[i], meta_files[i], dataset_names[i])
  }
  
  # Merge datasets
  merged_data <- merge_datasets(datasets)
  
  # Use relative abundance data (already converted)
  rel_abundance <- use_relative_abundance(merged_data$otu_data)
  
  # Identify differential microbes
  differential_results <- identify_differential_microbes(rel_abundance, merged_data$meta_data)
  
  # Create summary statistics
  summary_stats <- create_summary_stats(differential_results)
  
  # Create visualizations
  plots <- create_visualizations(differential_results, output_dir)
  
  # Save results
  cat("Saving results...\n")
  
  # Save differential microbes results
  write.csv(differential_results, file.path(output_dir, "differential_microbes.csv"), row.names = FALSE)
  
  # Save summary statistics
  summary_df <- data.frame(
    metric = names(summary_stats),
    value = unlist(summary_stats),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
  
  # Save detailed results
  sink(file.path(output_dir, "analysis_summary.txt"))
  cat("=== DIFFERENTIAL ABUNDANCE ANALYSIS SUMMARY ===\n\n")
  cat("Total datasets analyzed:", length(otu_files), "\n")
  cat("Dataset names:", paste(dataset_names, collapse = ", "), "\n")
  cat("Total samples:", ncol(merged_data$otu_data), "\n")
  cat("Total taxa before filtering:", nrow(merged_data$otu_data), "\n")
  cat("Total differential microbes identified:", summary_stats$total_differential_microbes, "\n")
  cat("Up-regulated microbes:", summary_stats$up_regulated, "\n")
  cat("Down-regulated microbes:", summary_stats$down_regulated, "\n")
  cat("Mean log2 fold change:", round(summary_stats$mean_log2_fc, 3), "\n")
  cat("Median log2 fold change:", round(summary_stats$median_log2_fc, 3), "\n")
  cat("Minimum q-value:", round(summary_stats$min_q_value, 6), "\n")
  cat("Maximum q-value:", round(summary_stats$max_q_value, 6), "\n")
  sink()
  
  # Print summary to console
  cat("\n=== SUMMARY ===\n")
  cat("Total differential microbes identified:", summary_stats$total_differential_microbes, "\n")
  cat("Up-regulated microbes:", summary_stats$up_regulated, "\n")
  cat("Down-regulated microbes:", summary_stats$down_regulated, "\n")
  cat("Mean log2 fold change:", round(summary_stats$mean_log2_fc, 3), "\n")
  
  cat("\nTop 10 differential microbes:\n")
  print(head(differential_results[, c("taxon", "log2_fold_change", "q_value")], 10))
  
  cat("\nAnalysis completed! Results saved to:", output_dir, "\n")
  
  return(list(
    differential_results = differential_results,
    summary_stats = summary_stats,
    plots = plots,
    merged_data = merged_data
  ))
}

# Command line interface
if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 6) {
    cat("Usage: Rscript differential_microbes.R <otu1> <meta1> <otu2> <meta2> <otu3> <meta3> [dataset1] [dataset2] [dataset3] [output_dir]\n")
    cat("Example: Rscript differential_microbes.R otu1.csv meta1.csv otu2.csv meta2.csv otu3.csv meta3.csv dataset1 dataset2 dataset3 results/\n")
    quit(status = 1)
  }
  
  otu_files <- args[1:3]
  meta_files <- args[4:6]
  dataset_names <- if (length(args) > 6) args[7:9] else paste0("dataset", 1:3)
  output_dir <- if (length(args) > 9) args[10] else "."
  
  # Check if files exist
  all_files <- c(otu_files, meta_files)
  for (file in all_files) {
    if (!file.exists(file)) {
      stop("Error: File not found:", file)
    }
  }
  
  # Run analysis
  results <- main(otu_files, meta_files, dataset_names, output_dir)
}

# Example usage (uncomment to use interactively):
# results <- main(
#   otu_files = c("otu1.csv", "otu2.csv", "otu3.csv"),
#   meta_files = c("meta1.csv", "meta2.csv", "meta3.csv"),
#   dataset_names = c("dataset1", "dataset2", "dataset3"),
#   output_dir = "differential_results"
# )
