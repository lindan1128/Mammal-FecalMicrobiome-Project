#!/usr/bin/env Rscript

# Differential abundance analysis script using NetMoss2
# Input: 3 OTU files and 3 metadata files
# Output: Differential microbes using NetMoss2

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(NetMoss)

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

# Function to calculate log2 fold change
calculate_log2_fold_change <- function(otu_data, meta_data) {
  cat("Calculating log2 fold change...\n")
  
  # Separate cases and controls
  case_samples <- rownames(meta_data)[meta_data$group == "sick"]
  control_samples <- rownames(meta_data)[meta_data$group == "normal"]
  
  cat("  Case samples:", length(case_samples), "\n")
  cat("  Control samples:", length(control_samples), "\n")
  
  # Calculate mean relative abundance for each group
  case_means <- rowMeans(otu_data[, case_samples, drop = FALSE], na.rm = TRUE)
  control_means <- rowMeans(otu_data[, control_samples, drop = FALSE], na.rm = TRUE)
  
  # Add small pseudocount to avoid log(0)
  pseudocount <- 1e-6
  case_means <- case_means + pseudocount
  control_means <- control_means + pseudocount
  
  # Calculate log2 fold change (case/control)
  log2_fc <- log2(case_means / control_means)
  
  return(log2_fc)
}

# Function to run NetMoss2 analysis
run_netmoss2_analysis <- function(otu_data, meta_data, output_dir, analysis_name) {
  cat("Running NetMoss2 analysis:", analysis_name, "\n")
  
  # Create integrated directory structure
  integrated_dir <- file.path(output_dir, paste0("integrated_", analysis_name, "_netmoss2"))
  if (dir.exists(integrated_dir)) {
    unlink(integrated_dir, recursive = TRUE)
  }
  dir.create(integrated_dir, recursive = TRUE)
  
  # Create subdirectories
  case_dir <- file.path(integrated_dir, "case_dir")
  control_dir <- file.path(integrated_dir, "control_dir")
  net_case_dir <- file.path(integrated_dir, "net_case_dir")
  net_control_dir <- file.path(integrated_dir, "net_control_dir")
  
  dir.create(case_dir, recursive = TRUE)
  dir.create(control_dir, recursive = TRUE)
  dir.create(net_case_dir, recursive = TRUE)
  dir.create(net_control_dir, recursive = TRUE)
  
  # Suppress warnings
  options(warn = -1)
  
  # Split data into case and control groups
  case_samples <- rownames(meta_data)[meta_data$group == "sick"]
  control_samples <- rownames(meta_data)[meta_data$group == "normal"]
  
  if (length(case_samples) == 0 || length(control_samples) == 0) {
    stop("No case or control samples found")
  }
  
  # Prepare case and control data
  case_data <- otu_data[, case_samples, drop = FALSE]
  control_data <- otu_data[, control_samples, drop = FALSE]
  
  cat("  Case data dimensions:", dim(case_data), "\n")
  cat("  Control data dimensions:", dim(control_data), "\n")
  
  # Prepare data for NetMoss2 format (taxa as rows, samples as columns)
  # Add empty first column name
  if (ncol(case_data) > 0) {
    colnames(case_data)[1] <- ""
  }
  if (ncol(control_data) > 0) {
    colnames(control_data)[1] <- ""
  }
  
  # Save to files
  case_file <- file.path(case_dir, "d_study1.txt")
  control_file <- file.path(control_dir, "h_study1.txt")
  
  write.table(case_data, case_file, sep = "\t", row.names = TRUE, quote = FALSE, col.names = TRUE, na = "")
  write.table(control_data, control_file, sep = "\t", row.names = TRUE, quote = FALSE, col.names = TRUE, na = "")
  
  # Build networks using Spearman correlation
  cat("  Building networks...\n")
  case_cor <- cor(t(case_data), method = "spearman")
  control_cor <- cor(t(control_data), method = "spearman")
  
  # Replace NA values with 0 and ensure diagonal is 1
  case_cor[is.na(case_cor)] <- 0
  control_cor[is.na(control_cor)] <- 0
  diag(case_cor) <- 1
  diag(control_cor) <- 1
  
  # Add empty first column name
  if (ncol(case_cor) > 0) {
    colnames(case_cor)[1] <- ""
  }
  if (ncol(control_cor) > 0) {
    colnames(control_cor)[1] <- ""
  }
  
  # Save network files
  net_case_file <- file.path(net_case_dir, "d_study_net1.txt")
  net_control_file <- file.path(net_control_dir, "h_study_net1.txt")
  
  write.table(case_cor, net_case_file, sep = "\t", row.names = TRUE, quote = FALSE, col.names = TRUE, na = "")
  write.table(control_cor, net_control_file, sep = "\t", row.names = TRUE, quote = FALSE, col.names = TRUE, na = "")
  
  # Run NetMoss2
  cat("  Running NetMoss2...\n")
  tryCatch({
    # Change to the integrated directory
    original_wd <- getwd()
    setwd(integrated_dir)
    
    # Read directory paths
    case_dir_path <- paste0(getwd(), "/case_dir")
    control_dir_path <- paste0(getwd(), "/control_dir")
    net_case_dir_path <- paste0(getwd(), "/net_case_dir")
    net_control_dir_path <- paste0(getwd(), "/net_control_dir")
    
    # Run NetMoss2
    integrated_results <- NetMoss(case_dir = case_dir_path,
                                 control_dir = control_dir_path,
                                 net_case_dir = net_case_dir_path,
                                 net_control_dir = net_control_dir_path)
    
    # Return to original directory
    setwd(original_wd)
    
    # Extract results
    netmoss_result <- integrated_results[[1]]
    
    # Create results dataframe
    netmoss_results <- data.frame(
      feature = netmoss_result$taxon_names,
      taxon_names = netmoss_result$taxon_names,
      NetMoss_Score = netmoss_result$NetMoss_Score,
      p.val = netmoss_result$p.val,
      p.adj = netmoss_result$p.adj,
      study_name = analysis_name,
      stringsAsFactors = FALSE
    )
    
    cat("  NetMoss2 analysis completed successfully\n")
    cat("  Total associations tested:", nrow(netmoss_results), "\n")
    
    return(netmoss_results)
    
  }, error = function(e) {
    # Return to original directory in case of error
    setwd(original_wd)
    cat("  Error in NetMoss2 analysis:", e$message, "\n")
    return(NULL)
  }, finally = {
    # Ensure we return to original directory
    setwd(original_wd)
    options(warn = 0)  # Restore warnings
  })
}

# Function to filter NetMoss2 results
filter_netmoss2_results <- function(netmoss_results, otu_data, meta_data, q_threshold = 0.05, netmoss_threshold = 0.3, log2fc_threshold = 1) {
  cat("Filtering NetMoss2 results...\n")
  cat("  Q-value threshold (p.adj):", q_threshold, "\n")
  cat("  NetMoss2 Score threshold:", netmoss_threshold, "\n")
  cat("  Log2 fold change threshold:", log2fc_threshold, "\n")
  
  if (is.null(netmoss_results)) {
    cat("  No results to filter\n")
    return(data.frame())
  }
  
  # Calculate log2 fold change
  log2_fc <- calculate_log2_fold_change(otu_data, meta_data)
  
  # Add log2 fold change to results
  netmoss_results$log2_fold_change <- log2_fc[netmoss_results$feature]
  
  # Apply selection criteria: p.adj < 0.05, NetMoss_Score > 0.3, and |log2_fold_change| > 1
  filtered_results <- netmoss_results %>% 
    filter(p.adj < q_threshold) %>%
    filter(NetMoss_Score > netmoss_threshold) %>%
    filter(abs(log2_fold_change) > log2fc_threshold)
  
  cat("  Results before filtering:", nrow(netmoss_results), "\n")
  cat("  Results after filtering:", nrow(filtered_results), "\n")
  
  # Add additional information
  if (nrow(filtered_results) > 0) {
    filtered_results$regulation <- ifelse(filtered_results$log2_fold_change > 0, "Up-regulated", "Down-regulated")
    filtered_results$abs_log2fc <- abs(filtered_results$log2_fold_change)
    
    # Sort by NetMoss2 Score
    filtered_results <- filtered_results[order(filtered_results$NetMoss_Score, decreasing = TRUE), ]
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
      mean_netmoss_score = NA,
      median_netmoss_score = NA,
      mean_log2_fc = NA,
      median_log2_fc = NA,
      min_q_value = NA,
      max_q_value = NA
    )
  } else {
    # Count up-regulated and down-regulated microbes
    up_regulated <- sum(filtered_results$log2_fold_change > 0)
    down_regulated <- sum(filtered_results$log2_fold_change < 0)
    
    # Summary statistics
    summary_stats <- list(
      analysis_name = analysis_name,
      total_differential_microbes = nrow(filtered_results),
      up_regulated = up_regulated,
      down_regulated = down_regulated,
      mean_netmoss_score = mean(filtered_results$NetMoss_Score),
      median_netmoss_score = median(filtered_results$NetMoss_Score),
      mean_log2_fc = mean(filtered_results$log2_fold_change),
      median_log2_fc = median(filtered_results$log2_fold_change),
      min_q_value = min(filtered_results$p.adj),
      max_q_value = max(filtered_results$p.adj)
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
  p1 <- ggplot(filtered_results, aes(x = log2_fold_change, y = -log10(p.adj))) +
    geom_point(aes(color = ifelse(p.adj < 0.05, "Significant", "Not significant")), alpha = 0.7) +
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
    
    p2 <- ggplot(top_microbes, aes(x = feature, y = NetMoss_Score, 
                                   fill = ifelse(log2_fold_change > 0, "Up-regulated", "Down-regulated"))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = paste("Top 20 Differential Microbes -", analysis_name),
           x = "Taxon",
           y = "NetMoss2 Score",
           fill = "Regulation") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8)) +
      scale_fill_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue"))
  } else {
    p2 <- NULL
  }
  
  # NetMoss2 Score vs Log2 Fold Change scatter plot
  p3 <- ggplot(filtered_results, aes(x = log2_fold_change, y = NetMoss_Score, 
                                     color = ifelse(p.adj < 0.05, "Significant", "Not significant"))) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "green") +
    labs(title = paste("NetMoss2 Score vs Log2 Fold Change -", analysis_name),
         x = "Log2 Fold Change (Case/Control)",
         y = "NetMoss2 Score",
         color = "Significance") +
    theme_minimal() +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray"))
  
  # Save plots
  ggsave(file.path(output_dir, paste0("volcano_plot_", analysis_name, ".png")), 
         p1, width = 10, height = 6, dpi = 300)
  
  if (!is.null(p2)) {
    ggsave(file.path(output_dir, paste0("top_differential_microbes_", analysis_name, ".png")), 
           p2, width = 10, height = 8, dpi = 300)
  }
  
  ggsave(file.path(output_dir, paste0("netmoss_score_vs_log2fc_", analysis_name, ".png")), 
         p3, width = 10, height = 6, dpi = 300)
  
  cat("  Plots saved to:", output_dir, "\n")
  
  return(list(volcano_plot = p1, top_microbes_plot = p2, netmoss_scatter_plot = p3))
}

# Main function
main <- function(otu_files, meta_files, dataset_names, output_dir = ".") {
  cat("Starting differential abundance analysis using NetMoss2...\n")
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
  
  # Run NetMoss2 analysis
  cat("\n=== NETMOSS2 ANALYSIS ===\n")
  netmoss_results <- run_netmoss2_analysis(filtered_otu_data, merged_data$meta_data, output_dir, "integrated")
  
  # Filter results
  filtered_results <- filter_netmoss2_results(netmoss_results, filtered_otu_data, merged_data$meta_data, 
                                             q_threshold = 0.05, netmoss_threshold = 0.3, log2fc_threshold = 1)
  
  # Create summary statistics
  summary_stats <- create_summary_stats(filtered_results, "NetMoss2")
  
  # Create visualizations
  plots <- create_visualizations(filtered_results, output_dir, "NetMoss2")
  
  # Save results
  cat("\nSaving results...\n")
  
  # Save filtered differential microbes results
  write.csv(filtered_results, file.path(output_dir, "differential_microbes_netmoss2.csv"), row.names = FALSE)
  
  # Save all NetMoss2 results
  if (!is.null(netmoss_results)) {
    write.csv(netmoss_results, file.path(output_dir, "netmoss2_all_results.csv"), row.names = FALSE)
  }
  
  # Save summary statistics
  summary_df <- data.frame(
    metric = names(summary_stats),
    value = unlist(summary_stats),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, file.path(output_dir, "summary_statistics_netmoss2.csv"), row.names = FALSE)
  
  # Save detailed results
  sink(file.path(output_dir, "netmoss2_analysis_summary.txt"))
  cat("=== NETMOSS2 DIFFERENTIAL ABUNDANCE ANALYSIS SUMMARY ===\n\n")
  cat("Total datasets analyzed:", length(otu_files), "\n")
  cat("Dataset names:", paste(dataset_names, collapse = ", "), "\n")
  cat("Total samples:", ncol(merged_data$otu_data), "\n")
  cat("Total taxa before filtering:", nrow(merged_data$otu_data), "\n")
  cat("Total taxa after filtering:", nrow(filtered_otu_data), "\n")
  cat("Total associations tested:", ifelse(is.null(netmoss_results), 0, nrow(netmoss_results)), "\n")
  cat("Total differential microbes identified:", summary_stats$total_differential_microbes, "\n")
  cat("Up-regulated microbes:", summary_stats$up_regulated, "\n")
  cat("Down-regulated microbes:", summary_stats$down_regulated, "\n")
  cat("Mean NetMoss2 Score:", round(summary_stats$mean_netmoss_score, 3), "\n")
  cat("Median NetMoss2 Score:", round(summary_stats$median_netmoss_score, 3), "\n")
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
  cat("Mean NetMoss2 Score:", round(summary_stats$mean_netmoss_score, 3), "\n")
  cat("Mean log2 fold change:", round(summary_stats$mean_log2_fc, 3), "\n")
  
  if (nrow(filtered_results) > 0) {
    cat("\nTop 10 differential microbes:\n")
    print(head(filtered_results[, c("feature", "NetMoss_Score", "log2_fold_change", "p.adj")], 10))
  } else {
    cat("\nNo differential microbes found with the specified criteria.\n")
  }
  
  cat("\nAnalysis completed! Results saved to:", output_dir, "\n")
  
  return(list(
    filtered_results = filtered_results,
    netmoss_results = netmoss_results,
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
    cat("Usage: Rscript netmoss2.R <otu1> <meta1> <otu2> <meta2> <otu3> <meta3> [dataset1] [dataset2] [dataset3] [output_dir]\n")
    cat("Example: Rscript netmoss2.R otu1.csv meta1.csv otu2.csv meta2.csv otu3.csv meta3.csv dataset1 dataset2 dataset3 results/\n")
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
#   output_dir = "netmoss2_results"
# )
