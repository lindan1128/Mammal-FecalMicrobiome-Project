# Required packages
library(MMUPHin)
library(Hmisc)
library(dplyr)
library(tidyr)
library(purrr)
library(circlize)
library(stringr)
library(readr)

otu_files <- list.files(pattern = 'otu_table_.*\\.tsv$') # List of OTU tables
meta_files <- list.files(pattern = 'metadata_.*\\.tsv$') # List of metadata files
output_prefix <- 'host_network'

# Read and batch-correct OTU tables
otu_list <- lapply(otu_files, read_tsv)
meta_list <- lapply(meta_files, read_tsv)

# Adjust batch effects for each study
otu_corrected_list <- map2(otu_list, meta_list, function(otu, meta) {
  otu <- as.data.frame(otu)
  rownames(otu) <- otu[[1]]
  otu <- otu[,-1]
  meta <- as.data.frame(meta)
  rownames(meta) <- meta[[1]]
  meta <- meta[,-1]
  common_samples <- intersect(colnames(otu), rownames(meta))
  otu <- otu[, common_samples]
  meta <- meta[common_samples, ]
  # Determine batch variable
  if ("project_id" %in% colnames(meta)) {
    batch_var <- "project_id"
  } else if ("breed" %in% colnames(meta)) {
    batch_var <- "breed"
  } else {
    stop("No suitable batch variable (project_id or breed) found in metadata.")
  }
  # Covariates: all except host, group, and batch
  covariate_vars <- setdiff(colnames(meta), c("host", "group", batch_var))
  if (length(covariate_vars) == 0) covariate_vars <- NULL
  adjusted <- adjust_batch(
    feature_abd = otu,
    batch = batch_var,
    covariates = covariate_vars,
    data = meta
  )
  as.data.frame(adjusted$feature_abd_adj)
})

# Merge all corrected OTU tables
all_otu <- reduce(otu_corrected_list, function(x, y) {
  full_join(x %>% rownames_to_column('taxa'), y %>% rownames_to_column('taxa'), by = 'taxa') %>%
    column_to_rownames('taxa')
})
all_otu[is.na(all_otu)] <- 0

# Merge all metadata
all_meta <- bind_rows(meta_list)
all_meta <- all_meta[match(colnames(all_otu), all_meta[[1]]), ]
rownames(all_meta) <- all_meta[[1]]
all_meta <- all_meta[,-1]

# Genus-level abundance (assume taxa names contain genus info, e.g., g__Genus)
genus_names <- sapply(strsplit(rownames(all_otu), "\\|"), function(x) {
  g <- x[grepl('^g__', x)]
  if(length(g) > 0) gsub('g__', '', g) else NA
})
all_otu_genus <- rowsum(as.matrix(all_otu), genus_names, na.rm=TRUE)

# For each host, calculate co-occurrence network
hosts <- unique(all_meta$host)
for (host in hosts) {
  host_samples <- rownames(all_meta)[all_meta$host == host & all_meta$group == 'normal']
  otu_host <- all_otu_genus[, host_samples, drop=FALSE]
  otu_host <- otu_host[rowSums(otu_host) > 0, ]
  otu_host <- t(otu_host)
  # Spearman correlation
  cor_res <- rcorr(as.matrix(otu_host), type = 'spearman')
  cor_mat <- cor_res$r
  p_mat <- cor_res$P
  # FDR correction
  p_mat[lower.tri(p_mat, diag=TRUE)] <- NA
  p_vec <- as.vector(p_mat)
  fdr_vec <- p.adjust(p_vec, method = 'fdr')
  fdr_mat <- matrix(NA, nrow=nrow(p_mat), ncol=ncol(p_mat))
  fdr_mat[upper.tri(fdr_mat)] <- fdr_vec
  # Filter significant correlations
  sig_idx <- which(abs(cor_mat) > 0.5 & fdr_mat < 0.05, arr.ind=TRUE)
  if(nrow(sig_idx) > 0) {
    links <- data.frame(
      from = rownames(cor_mat)[sig_idx[,1]],
      to = colnames(cor_mat)[sig_idx[,2]],
      rho = cor_mat[sig_idx],
      fdr = fdr_mat[sig_idx]
    )
    # Save network table
    write_tsv(links, paste0(output_prefix, '_', host, '_network.tsv'))
    # Plot network
    pdf(paste0(output_prefix, '_', host, '_chordDiagram.pdf'), width=8, height=8)
    chordDiagram(links[,c('from','to')], transparency = 0.5)
    dev.off()
  }
} 