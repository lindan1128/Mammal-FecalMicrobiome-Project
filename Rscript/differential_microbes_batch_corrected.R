#!/usr/bin/env Rscript

# Differential abundance analysis script with batch effect correction
# Input: 3 OTU files and 3 metadata files
# Output: Differential microbes after MMUPHin and sva ComBat batch correction

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(reshape2)
library(MMUPHin)
library(sva)

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
  
  # Check if Age and Sex columns exist
  if (!"Age" %in% colnames(meta_data)) {
    cat("  Warning: 'Age' column not found in metadata file:", meta_file, "\n")
    meta_data$Age <- NA
  }
  if (!"Sex" %in% colnames(meta_data)) {
    cat("  Warning: 'Sex' column not found in metadata file:", meta_file, "\n")
    meta_data$Sex <- NA
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
    
    # Add dataset and study columns to metadata
    dataset$meta_data$dataset <- dataset_name
    dataset$meta_data$study <- dataset_name  # For batch effect correction
    
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
  cat("Studies in merged data:", table(merged_meta$study), "\n")
  
  return(list(
    otu_data = merged_otu,
    meta_data = merged_meta
  ))
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

# Function to perform MMUPHin batch correction
perform_mmuphin_correction <- function(rel_abundance, meta_data) {
  cat("Performing MMUPHin batch correction...\n")
  
  # Transpose data for MMUPHin (samples as rows, features as columns)
  otu_relative_t <- t(rel_abundance)
  
  # Filter metadata to match samples
  meta_data_filtered <- meta_data[rownames(otu_relative_t), ]
  
  # Check for missing covariates
  if (all(is.na(meta_data_filtered$Age))) {
    cat("  Warning: All Age values are missing, using only Sex as covariate\n")
    covariates <- c('Sex')
  } else if (all(is.na(meta_data_filtered$Sex))) {
    cat("  Warning: All Sex values are missing, using only Age as covariate\n")
    covariates <- c('Age')
  } else if (all(is.na(meta_data_filtered$Age)) && all(is.na(meta_data_filtered$Sex))) {
    cat("  Warning: Both Age and Sex are missing, performing batch correction without covariates\n")
    covariates <- NULL
  } else {
    covariates <- c('Age', 'Sex')
  }
  
  # Perform MMUPHin batch correction
  tryCatch({
    ab_result <- adjust_batch(otu_relative_t, 
                             batch = 'study', 
                             data = meta_data_filtered,
                             covariates = covariates)
    
    # Transpose back to original format (features as rows, samples as columns)
    otu_relative_t_mm <- t(ab_result$feature_abd_adj)
    
    cat("  MMUPHin correction completed successfully\n")
    return(otu_relative_t_mm)
    
  }, error = function(e) {
    cat("  Error in MMUPHin correction:", e$message, "\n")
    cat("  Returning original data\n")
    return(rel_abundance)
  })
}

# Function to perform sva ComBat batch correction
perform_sva_correction <- function(rel_abundance, meta_data) {
  cat("Performing sva ComBat batch correction...\n")
  
  # Filter metadata to match samples
  meta_data_filtered <- meta_data[colnames(rel_abundance), ]
  
  # Create model matrix with group as the variable of interest
  mod <- model.matrix(~ group, data = meta_data_filtered)
  
  # Add covariates to model matrix if available
  if (!all(is.na(meta_data_filtered$Age)) && !all(is.na(meta_data_filtered$Sex))) {
    cat("  Including Age and Sex as covariates\n")
    mod <- model.matrix(~ group + Age + Sex, data = meta_data_filtered)
  } else if (!all(is.na(meta_data_filtered$Age))) {
    cat("  Including Age as covariate\n")
    mod <- model.matrix(~ group + Age, data = meta_data_filtered)
  } else if (!all(is.na(meta_data_filtered$Sex))) {
    cat("  Including Sex as covariate\n")
    mod <- model.matrix(~ group + Sex, data = meta_data_filtered)
  }
  
  # Perform ComBat batch correction
  tryCatch({
    otu_combat <- ComBat(dat = rel_abundance, 
                        batch = meta_data_filtered$study, 
                        mod = mod,
                        par.prior = TRUE,
                        prior.plots = FALSE)
    
    cat("  ComBat correction completed successfully\n")
    return(otu_combat)
    
  }, error = function(e) {
    cat("  Error in ComBat correction:", e$message, "\n")
    cat("  Returning original data\n")
    return(rel_abundance)
  })
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
identify_differential_microbes <- function(rel_abundance, meta_data, correction_method = "original") {
  cat("Identifying differential microbes (", correction_method, ")...\n")
  
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
    correction_method = correction_method,
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
create_summary_stats <- function(differential_results, correction_method) {
  cat("Creating summary statistics for", correction_method, "...\n")
  
  # Count up-regulated and down-regulated microbes
  up_regulated <- sum(differential_results$log2_fold_change > 0)
  down_regulated <- sum(differential_results$log2_fold_change < 0)
  
  # Summary statistics
  summary_stats <- list(
    correction_method = correction_method,
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
create_visualizations <- function(differential_results, output_dir, correction_method) {
  cat("Creating visualizations for", correction_method, "...\n")
  
  # Volcano plot
  p1 <- ggplot(differential_results, aes(x = log2_fold_change, y = -log10(q_value))) +
    geom_point(aes(color = ifelse(q_value < 0.05, "Significant", "Not significant")), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    labs(title = paste("Volcano Plot -", correction_method),
         x = "Log2 Fold Change (Case/Control)",
         y = "-Log10(Q-value)",
         color = "Significance") +
    theme_minimal() +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray"))
  
  # Top differential microbes bar plot
  top_microbes <- head(differential_results, 20)
  if (nrow(top_microbes) > 0) {
    top_microbes$taxon <- factor(top_microbes$taxon, levels = rev(top_microbes$taxon))
    
    p2 <- ggplot(top_microbes, aes(x = taxon, y = log2_fold_change, 
                                   fill = ifelse(log2_fold_change > 0, "Up-regulated", "Down-regulated"))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = paste("Top 20 Differential Microbes -", correction_method),
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
  ggsave(file.path(output_dir, paste0("volcano_plot_", correction_method, ".png")), 
         p1, width = 10, height = 6, dpi = 300)
  
  if (!is.null(p2)) {
    ggsave(file.path(output_dir, paste0("top_differential_microbes_", correction_method, ".png")), 
           p2, width = 10, height = 8, dpi = 300)
  }
  
  cat("  Plots saved to:", output_dir, "\n")
  
  return(list(volcano_plot = p1, top_microbes_plot = p2))
}

# Main function
main <- function(otu_files, meta_files, dataset_names, output_dir = ".") {
  cat("Starting differential abundance analysis with batch correction...\n")
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
  
  # Apply filtering criteria to original data
  filtered_rel_abundance <- apply_filtering_criteria(merged_data$otu_data, merged_data$meta_data)
  
  # Perform batch corrections
  cat("\n=== BATCH CORRECTION ===\n")
  
  # MMUPHin correction
  mmuphin_corrected <- perform_mmuphin_correction(filtered_rel_abundance, merged_data$meta_data)
  
  # sva ComBat correction
  sva_corrected <- perform_sva_correction(filtered_rel_abundance, merged_data$meta_data)
  
  # Identify differential microbes for each method
  cat("\n=== DIFFERENTIAL ANALYSIS ===\n")
  
  # Original data
  original_results <- identify_differential_microbes(filtered_rel_abundance, merged_data$meta_data, "original")
  
  # MMUPHin corrected
  mmuphin_results <- identify_differential_microbes(mmuphin_corrected, merged_data$meta_data, "MMUPHin")
  
  # sva ComBat corrected
  sva_results <- identify_differential_microbes(sva_corrected, merged_data$meta_data, "sva_ComBat")
  
  # Create summary statistics
  original_summary <- create_summary_stats(original_results, "original")
  mmuphin_summary <- create_summary_stats(mmuphin_results, "MMUPHin")
  sva_summary <- create_summary_stats(sva_results, "sva_ComBat")
  
  # Create visualizations
  original_plots <- create_visualizations(original_results, output_dir, "original")
  mmuphin_plots <- create_visualizations(mmuphin_results, output_dir, "MMUPHin")
  sva_plots <- create_visualizations(sva_results, output_dir, "sva_ComBat")
  
  # Combine all results
  all_results <- rbind(original_results, mmuphin_results, sva_results)
  all_summaries <- rbind(
    data.frame(original_summary, stringsAsFactors = FALSE),
    data.frame(mmuphin_summary, stringsAsFactors = FALSE),
    data.frame(sva_summary, stringsAsFactors = FALSE)
  )
  
  # Save results
  cat("\nSaving results...\n")
  
  # Save all differential microbes results
  write.csv(all_results, file.path(output_dir, "differential_microbes_all_methods.csv"), row.names = FALSE)
  
  # Save individual results
  write.csv(original_results, file.path(output_dir, "differential_microbes_original.csv"), row.names = FALSE)
  write.csv(mmuphin_results, file.path(output_dir, "differential_microbes_MMUPHin.csv"), row.names = FALSE)
  write.csv(sva_results, file.path(output_dir, "differential_microbes_sva_ComBat.csv"), row.names = FALSE)
  
  # Save summary statistics
  write.csv(all_summaries, file.path(output_dir, "summary_statistics_all_methods.csv"), row.names = FALSE)
  
  # Save detailed results
  sink(file.path(output_dir, "analysis_summary.txt"))
  cat("=== DIFFERENTIAL ABUNDANCE ANALYSIS WITH BATCH CORRECTION ===\n\n")
  cat("Total datasets analyzed:", length(otu_files), "\n")
  cat("Dataset names:", paste(dataset_names, collapse = ", "), "\n")
  cat("Total samples:", ncol(merged_data$otu_data), "\n")
  cat("Total taxa before filtering:", nrow(merged_data$otu_data), "\n\n")
  
  cat("=== RESULTS SUMMARY ===\n")
  for (i in 1:nrow(all_summaries)) {
    summary_row <- all_summaries[i, ]
    cat("Method:", summary_row$correction_method, "\n")
    cat("  Total differential microbes:", summary_row$total_differential_microbes, "\n")
    cat("  Up-regulated microbes:", summary_row$up_regulated, "\n")
    cat("  Down-regulated microbes:", summary_row$down_regulated, "\n")
    cat("  Mean log2 fold change:", round(summary_row$mean_log2_fc, 3), "\n")
    cat("  Median log2 fold change:", round(summary_row$median_log2_fc, 3), "\n")
    cat("  Min q-value:", round(summary_row$min_q_value, 6), "\n")
    cat("  Max q-value:", round(summary_row$max_q_value, 6), "\n\n")
  }
  sink()
  
  # Print summary to console
  cat("\n=== SUMMARY ===\n")
  for (i in 1:nrow(all_summaries)) {
    summary_row <- all_summaries[i, ]
    cat("Method:", summary_row$correction_method, "\n")
    cat("  Total differential microbes:", summary_row$total_differential_microbes, "\n")
    cat("  Up-regulated:", summary_row$up_regulated, "Down-regulated:", summary_row$down_regulated, "\n")
  }
  
  cat("\nAnalysis completed! Results saved to:", output_dir, "\n")
  
  return(list(
    all_results = all_results,
    all_summaries = all_summaries,
    original_results = original_results,
    mmuphin_results = mmuphin_results,
    sva_results = sva_results,
    merged_data = merged_data
  ))
}

# Command line interface
if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 6) {
    cat("Usage: Rscript differential_microbes_batch_corrected.R <otu1> <meta1> <otu2> <meta2> <otu3> <meta3> [dataset1] [dataset2] [dataset3] [output_dir]\n")
    cat("Example: Rscript differential_microbes_batch_corrected.R otu1.csv meta1.csv otu2.csv meta2.csv otu3.csv meta3.csv dataset1 dataset2 dataset3 results/\n")
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
#   output_dir = "differential_results_batch_corrected"
# )
