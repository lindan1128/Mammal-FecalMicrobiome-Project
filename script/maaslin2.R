#!/usr/bin/env Rscript

# Differential abundance analysis script using MaAsLin2
# Input: 3 OTU files and 3 metadata files
# Output: Differential microbes using MaAsLin2

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(Maaslin2)

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

# Function to apply filtering criteria
apply_filtering_criteria <- function(otu_data, meta_data) {
  cat("Applying filtering criteria...\n")
  
  # Criterion 1: Remove taxa with relative abundance < 0.001 across all samples
  cat("  Criterion 1: Filtering taxa with relative abundance < 0.001...\n")
  max_rel_abundance <- apply(otu_data, 1, max)
  taxa_to_keep <- max_rel_abundance >= 0.001
  
  filtered_otu_data <- otu_data[taxa_to_keep, , drop = FALSE]
  cat("  Taxa before filtering:", nrow(otu_data), "\n")
  cat("  Taxa after filtering:", nrow(filtered_otu_data), "\n")
  
  return(filtered_otu_data)
}

# Function to run MaAsLin2 analysis
run_maaslin2_analysis <- function(otu_data, meta_data, output_dir, analysis_name) {
  cat("Running MaAsLin2 analysis:", analysis_name, "\n")
  
  # Prepare data for MaAsLin2
  # MaAsLin2 expects samples as rows and features as columns
  otu_for_maaslin <- t(otu_data)
  
  # Ensure sample names match
  common_samples <- intersect(rownames(otu_for_maaslin), rownames(meta_data))
  if (length(common_samples) == 0) {
    stop("No common samples between OTU and metadata")
  }
  
  otu_for_maaslin <- otu_for_maaslin[common_samples, , drop = FALSE]
  meta_for_maaslin <- meta_data[common_samples, , drop = FALSE]
  
  cat("  Final data dimensions for MaAsLin2:", nrow(otu_for_maaslin), "samples x", ncol(otu_for_maaslin), "taxa\n")
  
  # Create output directory for this analysis
  maaslin_output_dir <- file.path(output_dir, paste0("maaslin2_", analysis_name))
  if (!dir.exists(maaslin_output_dir)) {
    dir.create(maaslin_output_dir, recursive = TRUE)
  }
  
  # Run MaAsLin2
  cat("  Running MaAsLin2...\n")
  tryCatch({
    maaslin2_result <- Maaslin2(
      input_data = otu_for_maaslin,
      input_metadata = meta_for_maaslin,
      output = maaslin_output_dir,
      fixed_effects = 'group',
      reference = c("group,normal"), 
      normalization = "NONE",
      transform = "NONE",
      analysis_method = "LM"
    )
    
    # Extract results
    maaslin2_results <- maaslin2_result$results
    
    cat("  MaAsLin2 analysis completed successfully\n")
    cat("  Total associations tested:", nrow(maaslin2_results), "\n")
    
    return(maaslin2_results)
    
  }, error = function(e) {
    cat("  Error in MaAsLin2 analysis:", e$message, "\n")
    return(NULL)
  })
}

# Function to filter MaAsLin2 results
filter_maaslin2_results <- function(maaslin2_results, q_threshold = 0.05, log2fc_threshold = 1) {
  cat("Filtering MaAsLin2 results...\n")
  cat("  Q-value threshold:", q_threshold, "\n")
  cat("  Log2 fold change threshold:", log2fc_threshold, "\n")
  
  if (is.null(maaslin2_results)) {
    cat("  No results to filter\n")
    return(data.frame())
  }
  
  # Apply selection criteria: Q-value < 0.05 and |log2 fold change| > 1
  filtered_results <- maaslin2_results %>% 
    filter(qval < q_threshold) %>%
    filter(abs(coef) > log2fc_threshold)
  
  cat("  Results before filtering:", nrow(maaslin2_results), "\n")
  cat("  Results after filtering:", nrow(filtered_results), "\n")
  
  # Add additional information
  if (nrow(filtered_results) > 0) {
    filtered_results$regulation <- ifelse(filtered_results$coef > 0, "Up-regulated", "Down-regulated")
    filtered_results$abs_log2fc <- abs(filtered_results$coef)
    
    # Sort by absolute log2 fold change
    filtered_results <- filtered_results[order(filtered_results$abs_log2fc, decreasing = TRUE), ]
  }
  
  return(filtered_results)
}

# Function to create summary statistics
create_summary_stats <- function(filtered_results, analysis_name) {
  cat("Creating summary statistics for", analysis_name, "...\n")
  
  if (nrow(filtered_results) == 0) {
    summary_stats <- list(
      analysis_name = analysis_name,
      total_differential_microbes = 0,
      up_regulated = 0,
      down_regulated = 0,
      mean_log2_fc = NA,
      median_log2_fc = NA,
      min_q_value = NA,
      max_q_value = NA
    )
  } else {
    # Count up-regulated and down-regulated microbes
    up_regulated <- sum(filtered_results$coef > 0)
    down_regulated <- sum(filtered_results$coef < 0)
    
    # Summary statistics
    summary_stats <- list(
      analysis_name = analysis_name,
      total_differential_microbes = nrow(filtered_results),
      up_regulated = up_regulated,
      down_regulated = down_regulated,
      mean_log2_fc = mean(filtered_results$coef),
      median_log2_fc = median(filtered_results$coef),
      min_q_value = min(filtered_results$qval),
      max_q_value = max(filtered_results$qval)
    )
  }
  
  return(summary_stats)
}

# Function to create visualizations
create_visualizations <- function(filtered_results, output_dir, analysis_name) {
  cat("Creating visualizations for", analysis_name, "...\n")
  
  if (nrow(filtered_results) == 0) {
    cat("  No results to visualize\n")
    return(list())
  }
  
  # Volcano plot
  p1 <- ggplot(filtered_results, aes(x = coef, y = -log10(qval))) +
    geom_point(aes(color = ifelse(qval < 0.05, "Significant", "Not significant")), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    labs(title = paste("Volcano Plot -", analysis_name),
         x = "Log2 Fold Change (Case/Control)",
         y = "-Log10(Q-value)",
         color = "Significance") +
    theme_minimal() +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray"))
  
  # Top differential microbes bar plot
  top_microbes <- head(filtered_results, 20)
  if (nrow(top_microbes) > 0) {
    top_microbes$feature <- factor(top_microbes$feature, levels = rev(top_microbes$feature))
    
    p2 <- ggplot(top_microbes, aes(x = feature, y = coef, 
                                   fill = ifelse(coef > 0, "Up-regulated", "Down-regulated"))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = paste("Top 20 Differential Microbes -", analysis_name),
           x = "Taxon",
           y = "Log2 Fold Change (Case/Control)",
           fill = "Regulation") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8)) +
      scale_fill_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue"))
  } else {
    p2 <- NULL
  }
  
  # Save plots
  ggsave(file.path(output_dir, paste0("volcano_plot_", analysis_name, ".png")), 
         p1, width = 10, height = 6, dpi = 300)
  
  if (!is.null(p2)) {
    ggsave(file.path(output_dir, paste0("top_differential_microbes_", analysis_name, ".png")), 
           p2, width = 10, height = 8, dpi = 300)
  }
  
  cat("  Plots saved to:", output_dir, "\n")
  
  return(list(volcano_plot = p1, top_microbes_plot = p2))
}

# Main function
main <- function(otu_files, meta_files, dataset_names, output_dir = ".") {
  cat("Starting differential abundance analysis using MaAsLin2...\n")
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
  
  # Apply filtering criteria
  filtered_otu_data <- apply_filtering_criteria(merged_data$otu_data, merged_data$meta_data)
  
  # Run MaAsLin2 analysis
  cat("\n=== MAASLIN2 ANALYSIS ===\n")
  maaslin2_results <- run_maaslin2_analysis(filtered_otu_data, merged_data$meta_data, output_dir, "integrated")
  
  # Filter results
  filtered_results <- filter_maaslin2_results(maaslin2_results, q_threshold = 0.05, log2fc_threshold = 1)
  
  # Create summary statistics
  summary_stats <- create_summary_stats(filtered_results, "MaAsLin2")
  
  # Create visualizations
  plots <- create_visualizations(filtered_results, output_dir, "MaAsLin2")
  
  # Save results
  cat("\nSaving results...\n")
  
  # Save filtered differential microbes results
  write.csv(filtered_results, file.path(output_dir, "differential_microbes_maaslin2.csv"), row.names = FALSE)
  
  # Save all MaAsLin2 results
  if (!is.null(maaslin2_results)) {
    write.csv(maaslin2_results, file.path(output_dir, "maaslin2_all_results.csv"), row.names = FALSE)
  }
  
  # Save summary statistics
  summary_df <- data.frame(
    metric = names(summary_stats),
    value = unlist(summary_stats),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, file.path(output_dir, "summary_statistics_maaslin2.csv"), row.names = FALSE)
  
  # Save detailed results
  sink(file.path(output_dir, "maaslin2_analysis_summary.txt"))
  cat("=== MAASLIN2 DIFFERENTIAL ABUNDANCE ANALYSIS SUMMARY ===\n\n")
  cat("Total datasets analyzed:", length(otu_files), "\n")
  cat("Dataset names:", paste(dataset_names, collapse = ", "), "\n")
  cat("Total samples:", ncol(merged_data$otu_data), "\n")
  cat("Total taxa before filtering:", nrow(merged_data$otu_data), "\n")
  cat("Total taxa after filtering:", nrow(filtered_otu_data), "\n")
  cat("Total associations tested:", ifelse(is.null(maaslin2_results), 0, nrow(maaslin2_results)), "\n")
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
  
  if (nrow(filtered_results) > 0) {
    cat("\nTop 10 differential microbes:\n")
    print(head(filtered_results[, c("feature", "coef", "qval")], 10))
  } else {
    cat("\nNo differential microbes found with the specified criteria.\n")
  }
  
  cat("\nAnalysis completed! Results saved to:", output_dir, "\n")
  
  return(list(
    filtered_results = filtered_results,
    maaslin2_results = maaslin2_results,
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
    cat("Usage: Rscript maaslin2.R <otu1> <meta1> <otu2> <meta2> <otu3> <meta3> [dataset1] [dataset2] [dataset3] [output_dir]\n")
    cat("Example: Rscript maaslin2.R otu1.csv meta1.csv otu2.csv meta2.csv otu3.csv meta3.csv dataset1 dataset2 dataset3 results/\n")
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
#   output_dir = "maaslin2_results"
# )
