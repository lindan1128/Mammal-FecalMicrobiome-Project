# Required packages
library(randomForest)
library(pROC)
library(dplyr)
library(caret)

otu_file <- 'otu_table.tsv'      # OTU table: rows are taxa, columns are samples
meta_file <- 'metadata.tsv'      # Metadata: rows are samples, columns include 'group'

# Read data
otu <- read.table(otu_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Ensure sample order consistency
common_samples <- intersect(colnames(otu), rownames(meta))
otu <- otu[, common_samples]
meta <- meta[common_samples, , drop = FALSE]

# Prepare data for modeling
X <- t(otu)
y <- as.factor(meta$group)
n_samples <- nrow(X)
sample_names <- rownames(X)

set.seed(123)
n_repeats <- 20
n_folds <- 3
auc_vec <- numeric(n_repeats)
gmhi_mat <- matrix(NA, nrow = n_samples, ncol = n_repeats * n_folds)
colnames(gmhi_mat) <- paste0('rep', rep(1:n_repeats, each = n_folds), '_fold', rep(1:n_folds, n_repeats))

for (rep in 1:n_repeats) {
  folds <- createFolds(y, k = n_folds, list = TRUE, returnTrain = FALSE)
  for (fold in 1:n_folds) {
    test_idx <- folds[[fold]]
    train_idx <- setdiff(1:n_samples, test_idx)
    rf_model <- randomForest(x = X[train_idx, , drop = FALSE], y = y[train_idx], ntree = 100)
    prob <- predict(rf_model, X[test_idx, , drop = FALSE], type = 'prob')[, 2]  # Probability for positive class
    gmhi_mat[test_idx, (rep - 1) * n_folds + fold] <- prob
    roc_obj <- roc(y[test_idx], prob, levels = rev(levels(y)))
    if (fold == 1) {
      auc_fold <- auc(roc_obj)
    } else {
      auc_fold <- c(auc_fold, auc(roc_obj))
    }
  }
  auc_vec[rep] <- mean(auc_fold)
}

mean_auc <- mean(auc_vec)
ci_auc <- quantile(auc_vec, c(0.025, 0.975))
mean_gmhi <- rowMeans(gmhi_mat, na.rm = TRUE)

result <- list(
  mean_auc = mean_auc,
  auc_95ci = ci_auc,
  mean_gmhi = data.frame(sample = sample_names, mean_gmhi = mean_gmhi)
)

print(result) 