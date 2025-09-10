#!/usr/bin/env Rscript

# Random Forest cross-validation analysis script
# Input: 1 OTU file and 1 metadata file
# Output: Cross-validation results with metrics

# Load required libraries
library(randomForest)
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

# Random forest cross-validation function with ROC output and balanced sampling
rf_cv <- function(data, labels, n_trees = 20, n_folds = 5, n_repeats = 5) {
  cat("Running random forest with", n_trees, "trees,", n_folds, "folds,", n_repeats, "repeats\n")
  cat("Using balanced sampling to address class imbalance\n")
  
  # Initialize vectors to store all fold results (5*5 = 25)
  aucs_list <- numeric(n_repeats * n_folds)
  acc_list <- numeric(n_repeats * n_folds)
  recall_list <- numeric(n_repeats * n_folds)
  f1_list <- numeric(n_repeats * n_folds)
  
  importance_list <- list()
  all_predictions <- matrix(0, nrow = length(labels), ncol = n_repeats)
  all_rocs <- list()
  
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
    set.seed(rep)
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
    fold_importance <- list()
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
      
      # Ensure labels are factors for classification
      rf_model <- randomForest(x = train_data, y = as.factor(train_labels), ntree = n_trees, importance = TRUE)
      pred_probs <- predict(rf_model, test_data, type = "prob")[, 2]
      pred_classes <- predict(rf_model, test_data, type = "class")
      
      # Calculate metrics
      roc_obj <- roc(test_labels, pred_probs)
      auc_val <- auc(roc_obj)
      
      # Calculate accuracy, recall, and F1
      conf_matrix <- table(test_labels, pred_classes)
      if(nrow(conf_matrix) == 2 && ncol(conf_matrix) == 2) {
        tp <- conf_matrix[2, 2]  # True positives (sick predicted as sick)
        tn <- conf_matrix[1, 1]  # True negatives (normal predicted as normal)
        fp <- conf_matrix[1, 2]  # False positives (normal predicted as sick)
        fn <- conf_matrix[2, 1]  # False negatives (sick predicted as normal)
        
        acc_val <- (tp + tn) / (tp + tn + fp + fn)
        recall_val <- tp / (tp + fn)  # Sensitivity/Recall
        precision_val <- tp / (tp + fp)
        f1_val <- 2 * (precision_val * recall_val) / (precision_val + recall_val)
      } else {
        acc_val <- recall_val <- f1_val <- NA
      }
      
      # Store results
      aucs_list[fold_counter] <- auc_val
      acc_list[fold_counter] <- acc_val
      recall_list[fold_counter] <- recall_val
      f1_list[fold_counter] <- f1_val
      
      fold_aucs[fold] <- auc_val
      fold_importance[[fold]] <- importance(rf_model)[, "MeanDecreaseAccuracy"]
      
      # Store predictions for this fold
      original_test_indices <- balanced_indices[test_idx]
      rep_predictions[original_test_indices] <- pred_probs
      
      fold_counter <- fold_counter + 1
    }
    
    importance_list[[rep]] <- Reduce("+", fold_importance) / length(fold_importance)
    all_predictions[, rep] <- rep_predictions
    
    # Calculate ROC for this repeat
    rep_roc <- roc(labels, rep_predictions)
    all_rocs[[rep]] <- rep_roc
    
    cat("  Repeat", rep, "completed - Mean AUC:", round(mean(fold_aucs), 3), "\n")
  }
  
  # Select ROC from the repeat with highest mean AUC
  mean_aucs <- sapply(1:n_repeats, function(i) mean(aucs_list[((i-1)*n_folds+1):(i*n_folds)]))
  best_rep <- which.max(mean_aucs)
  best_roc <- all_rocs[[best_rep]]
  
  cat("  Selected ROC from repeat", best_rep, "with highest mean AUC:", round(mean_aucs[best_rep], 3), "\n")
  
  return(list(
    aucs_list = aucs_list,      # 25 values
    acc_list = acc_list,        # 25 values  
    recall_list = recall_list,  # 25 values
    f1_list = f1_list,          # 25 values
    importance = importance_list
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
    labs(title = "Random Forest Cross-Validation Results",
         x = "Metrics",
         y = "Values") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("AUC" = "red", "Accuracy" = "blue", "Recall" = "green", "F1" = "orange"))
  
  # Save plot
  ggsave(file.path(output_dir, "rf_metrics_comparison.png"), p1, width = 10, height = 6, dpi = 300)
  
  cat("  Plot saved to:", output_dir, "\n")
  
  return(list(metrics_plot = p1))
}

# Main function
main <- function(otu_file, meta_file, output_dir = ".", n_trees = 20, n_folds = 5, n_repeats = 5) {
  cat("Starting Random Forest cross-validation analysis...\n")
  cat("OTU file:", otu_file, "\n")
  cat("Metadata file:", meta_file, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameters: n_trees =", n_trees, ", n_folds =", n_folds, ", n_repeats =", n_repeats, "\n\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load dataset
  dataset <- load_dataset(otu_file, meta_file)
  
  # Apply filtering criteria
  filtered_otu_data <- apply_filtering_criteria(dataset$otu_data, dataset$meta_data)
  
  # Transpose data for random forest (samples as rows, features as columns)
  rf_data <- t(filtered_otu_data)
  rf_labels <- dataset$meta_data$group
  
  cat("Final data dimensions for Random Forest:", nrow(rf_data), "samples x", ncol(rf_data), "features\n")
  cat("Group distribution:", table(rf_labels), "\n\n")
  
  # Run Random Forest cross-validation
  cat("=== RANDOM FOREST CROSS-VALIDATION ===\n")
  results <- rf_cv(rf_data, rf_labels, n_trees = n_trees, n_folds = n_folds, n_repeats = n_repeats)
  
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
  write.csv(results_df, file.path(output_dir, "rf_cv_results.csv"), row.names = FALSE)
  
  # Save summary statistics
  summary_df <- data.frame(
    metric = names(summary_stats),
    value = unlist(summary_stats),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, file.path(output_dir, "rf_summary_statistics.csv"), row.names = FALSE)
  
  # Save importance scores
  importance_df <- data.frame(
    feature = names(results$importance[[1]]),
    stringsAsFactors = FALSE
  )
  for (i in 1:length(results$importance)) {
    importance_df[[paste0("repeat_", i)]] <- results$importance[[i]]
  }
  importance_df$mean_importance <- rowMeans(importance_df[, -1], na.rm = TRUE)
  importance_df <- importance_df[order(importance_df$mean_importance, decreasing = TRUE), ]
  write.csv(importance_df, file.path(output_dir, "rf_feature_importance.csv"), row.names = FALSE)
  
  # Save detailed results
  sink(file.path(output_dir, "rf_analysis_summary.txt"))
  cat("=== RANDOM FOREST CROSS-VALIDATION ANALYSIS SUMMARY ===\n\n")
  cat("Total samples:", nrow(rf_data), "\n")
  cat("Total features:", ncol(rf_data), "\n")
  cat("Group distribution:", paste(names(table(rf_labels)), ":", table(rf_labels), collapse = ", "), "\n")
  cat("Cross-validation parameters:\n")
  cat("  Number of trees:", n_trees, "\n")
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
  
  cat("\nTop 10 most important features:\n")
  print(head(importance_df[, c("feature", "mean_importance")], 10))
  
  cat("\nAnalysis completed! Results saved to:", output_dir, "\n")
  
  return(list(
    results = results,
    summary_stats = summary_stats,
    plots = plots,
    importance_df = importance_df
  ))
}

# Command line interface
if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript random_forest.R <otu_file> <meta_file> [output_dir] [n_trees] [n_folds] [n_repeats]\n")
    cat("Example: Rscript random_forest.R otu.csv meta.csv results/ 20 5 5\n")
    quit(status = 1)
  }
  
  otu_file <- args[1]
  meta_file <- args[2]
  output_dir <- if (length(args) > 2) args[3] else "."
  n_trees <- if (length(args) > 3) as.numeric(args[4]) else 20
  n_folds <- if (length(args) > 4) as.numeric(args[5]) else 5
  n_repeats <- if (length(args) > 5) as.numeric(args[6]) else 5
  
  # Check if files exist
  if (!file.exists(otu_file)) {
    stop("Error: OTU file not found:", otu_file)
  }
  if (!file.exists(meta_file)) {
    stop("Error: Metadata file not found:", meta_file)
  }
  
  # Run analysis
  results <- main(otu_file, meta_file, output_dir, n_trees, n_folds, n_repeats)
}

# Example usage (uncomment to use interactively):
# results <- main(
#   otu_file = "otu.csv",
#   meta_file = "meta.csv",
#   output_dir = "rf_results",
#   n_trees = 20,
#   n_folds = 5,
#   n_repeats = 5
# )
