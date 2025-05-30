# Required packages
library(limma)
library(dplyr)
library(readr)
library(stringr)


otu_file <- 'otu_table.tsv'      # OTU table: rows are taxa, columns are samples
meta_file <- 'metadata.tsv'      # Metadata: rows are samples, columns are variables
output_file <- 'otu_table_corrected.tsv'
group_col <- 'group'             # Column name for group variable in metadata

# Read data
otu <- read_tsv(otu_file)
meta <- read_tsv(meta_file)

# Check row/column name matching
otu <- as.data.frame(otu)
rownames(otu) <- otu[[1]]
otu <- otu[,-1]
meta <- as.data.frame(meta)
rownames(meta) <- meta[[1]]
meta <- meta[,-1]

# Ensure sample order consistency
common_samples <- intersect(colnames(otu), rownames(meta))
otu <- otu[, common_samples]
meta <- meta[common_samples, ]

# Group information
group <- meta[[group_col]]

# Test significance of each variable between groups
test_results <- lapply(setdiff(colnames(meta), group_col), function(var) {
  x <- meta[[var]]
  if (is.character(x) || is.factor(x)) {
    # Fisher's exact test for qualitative variables
    tbl <- table(x, group)
    if (all(dim(tbl) > 1)) {
      p <- fisher.test(tbl)$p.value
    } else {
      p <- NA
    }
    type <- 'qualitative'
  } else {
    # Wilcoxon rank sum test for quantitative variables
    p <- tryCatch({
      wilcox.test(x ~ group)$p.value
    }, error = function(e) NA)
    type <- 'quantitative'
  }
  data.frame(var = var, p = p, type = type)
})
test_results <- bind_rows(test_results)

# Select significant variables
batch_vars <- test_results %>% filter(type == 'qualitative' & p < 0.05) %>% pull(var)
covariate_vars <- test_results %>% filter(type == 'quantitative' & p < 0.05) %>% pull(var)

# Construct batch and covariate matrices
batch <- if(length(batch_vars) > 0) meta[, batch_vars, drop=FALSE] else NULL
covariates <- if(length(covariate_vars) > 0) meta[, covariate_vars, drop=FALSE] else NULL

# Batch correction of OTU abundance
otu_mat <- as.matrix(otu)
design <- model.matrix(~ group)
corrected_otu <- removeBatchEffect(otu_mat, batch=batch, covariates=as.matrix(covariates), design=design)
corrected_otu <- as.data.frame(corrected_otu)

# Keep original row names
corrected_otu <- cbind(taxa=rownames(corrected_otu), corrected_otu)

# Output corrected OTU table
write_tsv(corrected_otu, output_file)

# Output significant variable list
write_tsv(test_results, 'covariate_test_results.tsv')

cat('Batch correction finished!\n') 