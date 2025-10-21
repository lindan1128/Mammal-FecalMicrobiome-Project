#!/usr/bin/env Rscript

# Differential abundance analysis script using ANCOM-BC2
# Input: 3 OTU files and 3 metadata files
# Output: Differential microbes using ANCOM-BC2

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(ANCOMBC)

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

# Function to apply filtering criteria (for absolute abundance data)
apply_filtering_criteria <- function(otu_data, meta_data) {
  cat("Applying filtering criteria for absolute abundance data...\n")
  
  # Criterion 1: Remove taxa with total count < 10 across all samples
  cat("  Criterion 1: Filtering taxa with total count < 10...\n")
  total_counts <- rowSums(otu_data)
  taxa_to_keep <- total_counts >= 10
  
  filtered_otu_data <- otu_data[taxa_to_keep, , drop = FALSE]
  cat("  Taxa before filtering:", nrow(otu_data), "\n")
  cat("  Taxa after filtering:", nrow(filtered_otu_data), "\n")
  
  return(filtered_otu_data)
}

# Function to run ANCOM-BC2 analysis
run_ancombc2_analysis <- function(otu_data, meta_data, output_dir, analysis_name) {
  cat("Running ANCOM-BC2 analysis:", analysis_name, "\n")
  
  # Prepare data for ANCOM-BC2
  # Ensure sample names match
  common_samples <- intersect(colnames(otu_data), rownames(meta_data))
  if (length(common_samples) == 0) {
    stop("No common samples between OTU and metadata")
  }
  
  otu_for_ancombc <- otu_data[, common_samples, drop = FALSE]
  meta_for_ancombc <- meta_data[common_samples, , drop = FALSE]
  
  cat("  Final data dimensions for ANCOM-BC2:", nrow(otu_for_ancombc), "taxa x", ncol(otu_for_ancombc), "samples\n")
  
  # Create phyloseq object for ANCOM-BC2
  cat("  Creating phyloseq object...\n")
  tryCatch({
    otu_table_obj <- otu_table(otu_for_ancombc, taxa_are_rows = TRUE)
    sample_data_obj <- sample_data(meta_for_ancombc)
    ps_obj <- phyloseq(otu_table_obj, sample_data_obj)
    
    cat("  Phyloseq object created successfully\n")
    
    # Run ANCOM-BC2
    cat("  Running ANCOM-BC2...\n")
    ancombc_result <- ancombc2(
      data = ps_obj,
      fix_formula = "group",
      rand_formula = NULL,
      p_adj_method = "BH",
      pseudo = 0,
      pseudo_sens = TRUE,
      prv_cut = 0.10,
      lib_cut = 1000,
      s0_perc = 0.90,
      group = "group",
      struc_zero = TRUE,
      neg_lb = TRUE,
      alpha = 0.05,
      n_cl = 2,
      verbose = TRUE,
      global = TRUE,
      pairwise = FALSE,
      dunnet = FALSE,
      trend = FALSE,
      iter_control = list(tolerance = 1e-2, max_iter = 20, verbose = FALSE),
      em_control = list(tolerance = 1e-5, max_iter = 20),
      lme_control = lme4::lmerControl(),
      mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
      trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)
    )
    
    # Process ANCOM-BC2 results
    cat("  Processing ANCOM-BC2 results...\n")
    
    # Extract results
    ancombc_results <- ancombc_result$res
    
    # Check if results exist and have the expected columns
    cat("  ANCOM-BC2 results dimensions:", nrow(ancombc_results), "x", ncol(ancombc_results), "\n")
    cat("  ANCOM-BC2 results columns:", paste(colnames(ancombc_results), collapse = ", "), "\n")
    
    # Add feature names as a column
    ancombc_results$feature <- rownames(ancombc_results)
    
    cat("  ANCOM-BC2 analysis completed successfully\n")
    cat("  Total associations tested:", nrow(ancombc_results), "\n")
    
    return(ancombc_results)
    
  }, error = function(e) {
    cat("  Error in ANCOM-BC2 analysis:", e$message, "\n")
    return(NULL)
  })
}

# Function to calculate log2 fold change
calculate_log2_fold_change <- function(otu_data, meta_data) {
  cat("Calculating log2 fold change...\n")
  
  # Separate cases and controls
  case_samples <- rownames(meta_data)[meta_data$group == "sick"]
  control_samples <- rownames(meta_data)[meta_data$group == "normal"]
  
  cat("  Case samples:", length(case_samples), "\n")
  cat("  Control samples:", length(control_samples), "\n")
  
  # Calculate mean abundance for each group
  case_means <- rowMeans(otu_data[, case_samples, drop = FALSE], na.rm = TRUE)
  control_means <- rowMeans(otu_data[, control_samples, drop = FALSE], na.rm = TRUE)
  
  # Add small pseudocount to avoid log(0)
  pseudocount <- 1
  case_means <- case_means + pseudocount
  control_means <- control_means + pseudocount
  
  # Calculate log2 fold change (case/control)
  log2_fc <- log2(case_means / control_means)
  
  return(log2_fc)
}

# Function to filter ANCOM-BC2 results
filter_ancombc2_results <- function(ancombc_results, otu_data, meta_data, q_threshold = 0.05, log2fc_threshold = 1) {
  cat("Filtering ANCOM-BC2 results...\n")
  cat("  Q-value threshold (p_groupsick):", q_threshold, "\n")
  cat("  Log2 fold change threshold:", log2fc_threshold, "\n")
  
  if (is.null(ancombc_results)) {
    cat("  No results to filter\n")
    return(data.frame())
  }
  
  # Calculate log2 fold change
  log2_fc <- calculate_log2_fold_change(otu_data, meta_data)
  
  # Add log2 fold change to results
  ancombc_results$log2_fold_change <- log2_fc[ancombc_results$feature]
  
  # Apply selection criteria: p_groupsick < 0.05 and |log2_fold_change| > 1
  filtered_results <- ancombc_results %>% 
    filter(p_groupsick < q_threshold) %>%
    filter(abs(log2_fold_change) > log2fc_threshold)
  
  cat("  Results before filtering:", nrow(ancombc_results), "\n")
  cat("  Results after filtering:", nrow(filtered_results), "\n")
  
  # Add additional information
  if (nrow(filtered_results) > 0) {
    filtered_results$regulation <- ifelse(filtered_results$log2_fold_change > 0, "Up-regulated", "Down-regulated")
    filtered_results$abs_log2fc <- abs(filtered_results$log2_fold_change)
    
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
    up_regulated <- sum(filtered_results$log2_fold_change > 0)
    down_regulated <- sum(filtered_results$log2_fold_change < 0)
    
    # Summary statistics
    summary_stats <- list(
      analysis_name = analysis_name,
      total_differential_microbes = nrow(filtered_results),
      up_regulated = up_regulated,
      down_regulated = down_regulated,
      mean_log2_fc = mean(filtered_results$log2_fold_change),
      median_log2_fc = median(filtered_results$log2_fold_change),
      min_q_value = min(filtered_results$p_groupsick),
      max_q_value = max(filtered_results$p_groupsick)
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
  p1 <- ggplot(filtered_results, aes(x = log2_fold_change, y = -log10(p_groupsick))) +
    geom_point(aes(color = ifelse(p_groupsick < 0.05, "Significant", "Not significant")), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    labs(title = paste("Volcano Plot -", analysis_name),
         x = "Log2 Fold Change (Case/Control)",
         y = "-Log10(P-value)",
         color = "Significance") +
    theme_minimal() +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray"))
  
  # Top differential microbes bar plot
  top_microbes <- head(filtered_results, 20)
  if (nrow(top_microbes) > 0) {
    top_microbes$feature <- factor(top_microbes$feature, levels = rev(top_microbes$feature))
    
    p2 <- ggplot(top_microbes, aes(x = feature, y = log2_fold_change, 
                                   fill = ifelse(log2_fold_change > 0, "Up-regulated", "Down-regulated"))) +
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
  
  # Log2 fold change distribution plot
  p3 <- ggplot(filtered_results, aes(x = log2_fold_change, fill = regulation)) +
    geom_histogram(bins = 30, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    labs(title = paste("Log2 Fold Change Distribution -", analysis_name),
         x = "Log2 Fold Change (Case/Control)",
         y = "Count",
         fill = "Regulation") +
    theme_minimal() +
    scale_fill_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue"))
  
  # Save plots
  ggsave(file.path(output_dir, paste0("volcano_plot_", analysis_name, ".png")), 
         p1, width = 10, height = 6, dpi = 300)
  
  if (!is.null(p2)) {
    ggsave(file.path(output_dir, paste0("top_differential_microbes_", analysis_name, ".png")), 
           p2, width = 10, height = 8, dpi = 300)
  }
  
  ggsave(file.path(output_dir, paste0("log2fc_distribution_", analysis_name, ".png")), 
         p3, width = 10, height = 6, dpi = 300)
  
  cat("  Plots saved to:", output_dir, "\n")
  
  return(list(volcano_plot = p1, top_microbes_plot = p2, log2fc_distribution_plot = p3))
}

# Main function
main <- function(otu_files, meta_files, dataset_names, output_dir = ".") {
  cat("Starting differential abundance analysis using ANCOM-BC2...\n")
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
  
  # Run ANCOM-BC2 analysis
  cat("\n=== ANCOM-BC2 ANALYSIS ===\n")
  ancombc_results <- run_ancombc2_analysis(filtered_otu_data, merged_data$meta_data, output_dir, "integrated")
  
  # Filter results
  filtered_results <- filter_ancombc2_results(ancombc_results, filtered_otu_data, merged_data$meta_data, q_threshold = 0.05, log2fc_threshold = 1)
  
  # Create summary statistics
  summary_stats <- create_summary_stats(filtered_results, "ANCOM-BC2")
  
  # Create visualizations
  plots <- create_visualizations(filtered_results, output_dir, "ANCOM-BC2")
  
  # Save results
  cat("\nSaving results...\n")
  
  # Save filtered differential microbes results
  write.csv(filtered_results, file.path(output_dir, "differential_microbes_ancombc2.csv"), row.names = FALSE)
  
  # Save all ANCOM-BC2 results
  if (!is.null(ancombc_results)) {
    write.csv(ancombc_results, file.path(output_dir, "ancombc2_all_results.csv"), row.names = FALSE)
  }
  
  # Save summary statistics
  summary_df <- data.frame(
    metric = names(summary_stats),
    value = unlist(summary_stats),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, file.path(output_dir, "summary_statistics_ancombc2.csv"), row.names = FALSE)
  
  # Save detailed results
  sink(file.path(output_dir, "ancombc2_analysis_summary.txt"))
  cat("=== ANCOM-BC2 DIFFERENTIAL ABUNDANCE ANALYSIS SUMMARY ===\n\n")
  cat("Total datasets analyzed:", length(otu_files), "\n")
  cat("Dataset names:", paste(dataset_names, collapse = ", "), "\n")
  cat("Total samples:", ncol(merged_data$otu_data), "\n")
  cat("Total taxa before filtering:", nrow(merged_data$otu_data), "\n")
  cat("Total taxa after filtering:", nrow(filtered_otu_data), "\n")
  cat("Total associations tested:", ifelse(is.null(ancombc_results), 0, nrow(ancombc_results)), "\n")
  cat("Total differential microbes identified:", summary_stats$total_differential_microbes, "\n")
  cat("Up-regulated microbes:", summary_stats$up_regulated, "\n")
  cat("Down-regulated microbes:", summary_stats$down_regulated, "\n")
  cat("Mean log2 fold change:", round(summary_stats$mean_log2_fc, 3), "\n")
  cat("Median log2 fold change:", round(summary_stats$median_log2_fc, 3), "\n")
  cat("Minimum p-value:", round(summary_stats$min_q_value, 6), "\n")
  cat("Maximum p-value:", round(summary_stats$max_q_value, 6), "\n")
  sink()
  
  # Print summary to console
  cat("\n=== SUMMARY ===\n")
  cat("Total differential microbes identified:", summary_stats$total_differential_microbes, "\n")
  cat("Up-regulated microbes:", summary_stats$up_regulated, "\n")
  cat("Down-regulated microbes:", summary_stats$down_regulated, "\n")
  cat("Mean log2 fold change:", round(summary_stats$mean_log2_fc, 3), "\n")
  
  if (nrow(filtered_results) > 0) {
    cat("\nTop 10 differential microbes:\n")
    print(head(filtered_results[, c("feature", "log2_fold_change", "p_groupsick")], 10))
  } else {
    cat("\nNo differential microbes found with the specified criteria.\n")
  }
  
  cat("\nAnalysis completed! Results saved to:", output_dir, "\n")
  
  return(list(
    filtered_results = filtered_results,
    ancombc_results = ancombc_results,
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
    cat("Usage: Rscript ancombc2.R <otu1> <meta1> <otu2> <meta2> <otu3> <meta3> [dataset1] [dataset2] [dataset3] [output_dir]\n")
    cat("Example: Rscript ancombc2.R otu1.csv meta1.csv otu2.csv meta2.csv otu3.csv meta3.csv dataset1 dataset2 dataset3 results/\n")
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
#   output_dir = "ancombc2_results"
# )
