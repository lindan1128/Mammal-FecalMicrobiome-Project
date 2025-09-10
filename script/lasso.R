#!/usr/bin/env Rscript

# Lasso regression cross-validation analysis script
# Input: 1 OTU file and 1 metadata file
# Output: Cross-validation results with metrics

# Load required libraries
library(glmnet)
library(caret)
library(pROC)
library(dplyr)
library(ggplot2)

# Function to load and process dataset
load_dataset <- function(otu_file, meta_file) {
  cat("Loading dataset...\n")
  
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
    meta_data = meta_data
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

# Lasso cross-validation function with balanced sampling and best ROC selection
lasso_cv <- function(data, labels, n_folds = 5, n_repeats = 5) {
  cat("Running Lasso regression with", n_folds, "folds,", n_repeats, "repeats\n")
  cat("Using balanced sampling to address class imbalance\n")
  
  # Initialize vectors to store all fold results (5*5 = 25)
  aucs_list <- numeric(n_repeats * n_folds)
  acc_list <- numeric(n_repeats * n_folds)
  recall_list <- numeric(n_repeats * n_folds)
  f1_list <- numeric(n_repeats * n_folds)
  
  # Also keep the original summary metrics
  aucs <- numeric(n_repeats)
  all_predictions <- matrix(0, nrow = length(labels), ncol = n_repeats)
  all_rocs <- list()  # Store ROC objects for each repeat
  
  fold_counter <- 1  # Counter for all folds across repeats
  
  for(rep in 1:n_repeats) {
    set.seed(rep)
    
    # Balance the dataset before creating folds
    cat("  Repeat", rep, "- Balancing dataset...\n")
    
    # Get indices for each group
    normal_indices <- which(labels == "normal")
    sick_indices <- which(labels == "sick")
    
    # Determine which group is larger and needs downsampling
    if(length(normal_indices) > length(sick_indices)) {
      # Normal group is larger, downsample it
      major_indices <- normal_indices
      minor_indices <- sick_indices
      major_group <- "normal"
      minor_group <- "sick"
    } else {
      # Sick group is larger, downsample it
      major_indices <- sick_indices
      minor_indices <- normal_indices
      major_group <- "sick"
      minor_group <- "normal"
    }
    
    # Randomly sample from major group to match minor group size
    set.seed(rep)  # Ensure reproducible sampling
    balanced_major_indices <- sample(major_indices, size = length(minor_indices), replace = FALSE)
    
    # Create balanced dataset
    balanced_indices <- c(balanced_major_indices, minor_indices)
    balanced_data <- data[balanced_indices, ]
    balanced_labels <- labels[balanced_indices]
    
    cat("    Original - Normal:", length(normal_indices), "Sick:", length(sick_indices), "\n")
    cat("    Balanced - Normal:", sum(balanced_labels == "normal"), "Sick:", sum(balanced_labels == "sick"), "\n")
    
    # Create folds on balanced data
    folds <- createFolds(balanced_labels, k = n_folds)
    fold_aucs <- numeric(n_folds)
    rep_predictions <- numeric(length(labels))
    
    for(fold in 1:n_folds) {
      train_idx <- unlist(folds[-fold])
      test_idx <- unlist(folds[fold])
      
      train_data <- balanced_data[train_idx, ]
      train_labels <- balanced_labels[train_idx]
      test_data <- balanced_data[test_idx, ]
      test_labels <- balanced_labels[test_idx]
      
      cat("    Fold", fold, "- Training samples:", length(train_idx), "Test samples:", length(test_idx), "\n")
      cat("      Training - Normal:", sum(train_labels == "normal"), "Sick:", sum(train_labels == "sick"), "\n")
      cat("      Test - Normal:", sum(test_labels == "normal"), "Sick:", sum(test_labels == "sick"), "\n")
      
      # Convert labels to numeric for Lasso (0 = normal, 1 = sick)
      train_labels_numeric <- as.numeric(train_labels == "sick")
      test_labels_numeric <- as.numeric(test_labels == "sick")
      
      # Lasso training with adaptive strategy based on training set size
      if(nrow(train_data) < 10) {
        # Training set too small, use fixed lambda
        cat("      Warning: Training set too small (", nrow(train_data), "), using fixed lambda\n")
        fit <- glmnet(as.matrix(train_data), train_labels_numeric, 
                      family = "binomial", alpha = 1, lambda = 0.01)
        pred_probs <- predict(fit, as.matrix(test_data), s = 0.01, type = "response")
      } else {
        # Normal case: use cross-validation, ensure nfolds >= 4
        nfolds_to_use <- max(4, min(10, nrow(train_data) %/% 2))
        cat("      Using", nfolds_to_use, "-fold CV\n")
        cv_fit <- cv.glmnet(as.matrix(train_data), train_labels_numeric, 
                            family = "binomial", alpha = 1, nfolds = nfolds_to_use)
        pred_probs <- predict(cv_fit, as.matrix(test_data), s = "lambda.min", type = "response")
      }
      
      pred_classes <- ifelse(pred_probs > 0.5, 1, 0)
      
      # Calculate AUC
      roc_obj <- roc(test_labels_numeric, as.numeric(pred_probs))
      fold_aucs[fold] <- auc(roc_obj)
      
      # Confusion matrix calculation - force creation of 2x2 matrix
      # Ensure all classes appear in confusion matrix
      all_classes <- c(0, 1)
      
      # Create complete confusion matrix
      conf_matrix <- table(factor(test_labels_numeric, levels = all_classes), 
                          factor(pred_classes, levels = all_classes))
      
      cat("      Confusion Matrix:\n")
      print(conf_matrix)
      
      # Now confusion matrix is always 2x2
      # 0 at row 1 col 1, 1 at row 2 col 2
      tn <- conf_matrix[1, 1]  # True negatives (0 predicted as 0)
      fp <- conf_matrix[1, 2]  # False positives (0 predicted as 1)
      fn <- conf_matrix[2, 1]  # False negatives (1 predicted as 0)
      tp <- conf_matrix[2, 2]  # True positives (1 predicted as 1)
      
      acc_val <- (tp + tn) / (tp + tn + fp + fn)
      recall_val <- tp / (tp + fn)  # Sensitivity/Recall
      precision_val <- tp / (tp + fp)
      f1_val <- 2 * (precision_val * recall_val) / (precision_val + recall_val)
      
      cat("      TP:", tp, "TN:", tn, "FP:", fp, "FN:", fn, "\n")
      cat("      Precision:", round(precision_val, 3), "Recall:", round(recall_val, 3), "F1:", round(f1_val, 3), "\n")
      
      # Handle NaN values
      if(is.nan(f1_val)) f1_val <- NA
      if(is.nan(acc_val)) acc_val <- NA
      if(is.nan(recall_val)) recall_val <- NA
      
      # Add debugging information
      cat("      Test labels:", table(test_labels_numeric), "\n")
      cat("      Pred classes:", table(pred_classes), "\n")
      cat("      Pred probs range:", round(range(pred_probs), 3), "\n")
      cat("      Pred probs mean:", round(mean(pred_probs), 3), "\n")
      
      # Store results for this fold
      aucs_list[fold_counter] <- fold_aucs[fold]
      acc_list[fold_counter] <- acc_val
      recall_list[fold_counter] <- recall_val
      f1_list[fold_counter] <- f1_val
      
      # Store predictions for this fold (map back to original sample indices)
      original_test_indices <- balanced_indices[test_idx]
      rep_predictions[original_test_indices] <- as.numeric(pred_probs)
      
      cat("      Final - AUC:", round(fold_aucs[fold], 3), "Accuracy:", round(acc_val, 3), 
          "Recall:", round(recall_val, 3), "F1:", round(f1_val, 3), "\n")
      
      fold_counter <- fold_counter + 1
    }
    
    aucs[rep] <- mean(fold_aucs)
    all_predictions[, rep] <- rep_predictions
    
    # Calculate ROC for this repeat
    rep_roc <- roc(as.numeric(labels == "sick"), rep_predictions)
    all_rocs[[rep]] <- rep_roc
    
    cat("  Repeat", rep, "completed - Mean AUC:", round(aucs[rep], 3), "\n")
  }
  
  # Select ROC from the repeat with highest AUC
  best_rep <- which.max(aucs)
  best_roc <- all_rocs[[best_rep]]
  
  cat("  Selected ROC from repeat", best_rep, "with highest AUC:", round(aucs[best_rep], 3), "\n")
  
  return(list(
    aucs_list = aucs_list,      # 25 values
    acc_list = acc_list,        # 25 values  
    recall_list = recall_list,  # 25 values
    f1_list = f1_list,          # 25 values
    importance = list()  # Lasso doesn't have feature importance like RF
  ))
}

# Function to create summary statistics
create_summary_stats <- function(results) {
  cat("Creating summary statistics...\n")
  
  # Calculate summary statistics
  summary_stats <- list(
    mean_auc = mean(results$aucs_list, na.rm = TRUE),
    sd_auc = sd(results$aucs_list, na.rm = TRUE),
    mean_accuracy = mean(results$acc_list, na.rm = TRUE),
    sd_accuracy = sd(results$acc_list, na.rm = TRUE),
    mean_recall = mean(results$recall_list, na.rm = TRUE),
    sd_recall = sd(results$recall_list, na.rm = TRUE),
    mean_f1 = mean(results$f1_list, na.rm = TRUE),
    sd_f1 = sd(results$f1_list, na.rm = TRUE)
  )
  
  return(summary_stats)
}

# Function to create visualizations
create_visualizations <- function(results, output_dir) {
  cat("Creating visualizations...\n")
  
  # Create metrics comparison plot
  metrics_data <- data.frame(
    Metric = rep(c("AUC", "Accuracy", "Recall", "F1"), each = 25),
    Value = c(results$aucs_list, results$acc_list, results$recall_list, results$f1_list)
  )
  
  p1 <- ggplot(metrics_data, aes(x = Metric, y = Value, fill = Metric)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = "Lasso Cross-Validation Results",
         x = "Metrics",
         y = "Values") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("AUC" = "red", "Accuracy" = "blue", "Recall" = "green", "F1" = "orange"))
  
  # Save plot
  ggsave(file.path(output_dir, "lasso_metrics_comparison.png"), p1, width = 10, height = 6, dpi = 300)
  
  cat("  Plot saved to:", output_dir, "\n")
  
  return(list(metrics_plot = p1))
}

# Main function
main <- function(otu_file, meta_file, output_dir = ".", n_folds = 5, n_repeats = 5) {
  cat("Starting Lasso regression cross-validation analysis...\n")
  cat("OTU file:", otu_file, "\n")
  cat("Metadata file:", meta_file, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameters: n_folds =", n_folds, ", n_repeats =", n_repeats, "\n\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load dataset
  dataset <- load_dataset(otu_file, meta_file)
  
  # Apply filtering criteria
  filtered_otu_data <- apply_filtering_criteria(dataset$otu_data, dataset$meta_data)
  
  # Transpose data for lasso (samples as rows, features as columns)
  lasso_data <- t(filtered_otu_data)
  lasso_labels <- dataset$meta_data$group
  
  cat("Final data dimensions for Lasso:", nrow(lasso_data), "samples x", ncol(lasso_data), "features\n")
  cat("Group distribution:", table(lasso_labels), "\n\n")
  
  # Run Lasso cross-validation
  cat("=== LASSO CROSS-VALIDATION ===\n")
  results <- lasso_cv(lasso_data, lasso_labels, n_folds = n_folds, n_repeats = n_repeats)
  
  # Create summary statistics
  summary_stats <- create_summary_stats(results)
  
  # Create visualizations
  plots <- create_visualizations(results, output_dir)
  
  # Save results
  cat("\nSaving results...\n")
  
  # Save detailed results
  results_df <- data.frame(
    fold = rep(1:n_folds, n_repeats),
    repeat_num = rep(1:n_repeats, each = n_folds),
    auc = results$aucs_list,
    accuracy = results$acc_list,
    recall = results$recall_list,
    f1 = results$f1_list,
    stringsAsFactors = FALSE
  )
  write.csv(results_df, file.path(output_dir, "lasso_cv_results.csv"), row.names = FALSE)
  
  # Save summary statistics
  summary_df <- data.frame(
    metric = names(summary_stats),
    value = unlist(summary_stats),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, file.path(output_dir, "lasso_summary_statistics.csv"), row.names = FALSE)
  
  # Save detailed results
  sink(file.path(output_dir, "lasso_analysis_summary.txt"))
  cat("=== LASSO CROSS-VALIDATION ANALYSIS SUMMARY ===\n\n")
  cat("Total samples:", nrow(lasso_data), "\n")
  cat("Total features:", ncol(lasso_data), "\n")
  cat("Group distribution:", paste(names(table(lasso_labels)), ":", table(lasso_labels), collapse = ", "), "\n")
  cat("Cross-validation parameters:\n")
  cat("  Number of folds:", n_folds, "\n")
  cat("  Number of repeats:", n_repeats, "\n")
  cat("  Total iterations:", n_folds * n_repeats, "\n\n")
  cat("Results summary:\n")
  cat("  Mean AUC:", round(summary_stats$mean_auc, 3), "±", round(summary_stats$sd_auc, 3), "\n")
  cat("  Mean Accuracy:", round(summary_stats$mean_accuracy, 3), "±", round(summary_stats$sd_accuracy, 3), "\n")
  cat("  Mean Recall:", round(summary_stats$mean_recall, 3), "±", round(summary_stats$sd_recall, 3), "\n")
  cat("  Mean F1:", round(summary_stats$mean_f1, 3), "±", round(summary_stats$sd_f1, 3), "\n")
  sink()
  
  # Print summary to console
  cat("\n=== SUMMARY ===\n")
  cat("Mean AUC:", round(summary_stats$mean_auc, 3), "±", round(summary_stats$sd_auc, 3), "\n")
  cat("Mean Accuracy:", round(summary_stats$mean_accuracy, 3), "±", round(summary_stats$sd_accuracy, 3), "\n")
  cat("Mean Recall:", round(summary_stats$mean_recall, 3), "±", round(summary_stats$sd_recall, 3), "\n")
  cat("Mean F1:", round(summary_stats$mean_f1, 3), "±", round(summary_stats$sd_f1, 3), "\n")
  
  cat("\nAnalysis completed! Results saved to:", output_dir, "\n")
  
  return(list(
    results = results,
    summary_stats = summary_stats,
    plots = plots
  ))
}

# Command line interface
if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript lasso.R <otu_file> <meta_file> [output_dir] [n_folds] [n_repeats]\n")
    cat("Example: Rscript lasso.R otu.csv meta.csv results/ 5 5\n")
    quit(status = 1)
  }
  
  otu_file <- args[1]
  meta_file <- args[2]
  output_dir <- if (length(args) > 2) args[3] else "."
  n_folds <- if (length(args) > 3) as.numeric(args[4]) else 5
  n_repeats <- if (length(args) > 4) as.numeric(args[5]) else 5
  
  # Check if files exist
  if (!file.exists(otu_file)) {
    stop("Error: OTU file not found:", otu_file)
  }
  if (!file.exists(meta_file)) {
    stop("Error: Metadata file not found:", meta_file)
  }
  
  # Run analysis
  results <- main(otu_file, meta_file, output_dir, n_folds, n_repeats)
}

# Example usage (uncomment to use interactively):
# results <- main(
#   otu_file = "otu.csv",
#   meta_file = "meta.csv",
#   output_dir = "lasso_results",
#   n_folds = 5,
#   n_repeats = 5
# )
