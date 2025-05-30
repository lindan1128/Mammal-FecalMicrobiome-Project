# Required packages
library(Hmisc)
library(igraph)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

otu_file <- 'otu_table.tsv'      # OTU table: rows are ASVs, columns are samples
meta_file <- 'metadata.tsv'      # Metadata: rows are samples, columns include 'group'
output_file <- 'bootstrap_network_results.rds'

# Read data
otu <- read_tsv(otu_file)
meta <- read_tsv(meta_file)

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

# Function to perform bootstrapped correlation analysis for a group
group_bootstrap_correlation <- function(otu_group, n_boot = 100, prop = 0.9, rho_cut = 0.3, fdr_cut = 0.05) {
  n_samples <- ncol(otu_group)
  n_select <- ceiling(n_samples * prop)
  asv_names <- rownames(otu_group)
  cor_list <- vector("list", n_boot)
  for (i in seq_len(n_boot)) {
    set.seed(i)
    sel_samples <- sample(colnames(otu_group), n_select, replace = FALSE)
    otu_sub <- otu_group[, sel_samples, drop=FALSE]
    otu_sub <- otu_sub[rowSums(otu_sub) > 0, ]
    otu_sub <- t(otu_sub)
    cor_res <- rcorr(as.matrix(otu_sub), type = 'spearman')
    cor_mat <- cor_res$r
    p_mat <- cor_res$P
    p_mat[lower.tri(p_mat, diag=TRUE)] <- NA
    p_vec <- as.vector(p_mat)
    fdr_vec <- p.adjust(p_vec, method = 'fdr')
    fdr_mat <- matrix(NA, nrow=nrow(p_mat), ncol=ncol(p_mat))
    fdr_mat[upper.tri(fdr_mat)] <- fdr_vec
    sig_idx <- which(abs(cor_mat) > rho_cut & fdr_mat < fdr_cut, arr.ind=TRUE)
    if(nrow(sig_idx) > 0) {
      links <- data.frame(
        from = rownames(cor_mat)[sig_idx[,1]],
        to = colnames(cor_mat)[sig_idx[,2]],
        rho = cor_mat[sig_idx],
        fdr = fdr_mat[sig_idx]
      )
      cor_list[[i]] <- links
    } else {
      cor_list[[i]] <- data.frame(from=character(), to=character(), rho=numeric(), fdr=numeric())
    }
  }
  cor_list
}

# Function to get bootstrap network edges (present in at least 50% of bootstraps)
bootstrap_edges <- function(cor_list, min_boot = 0.5) {
  all_edges <- bind_rows(cor_list, .id = "boot")
  edge_counts <- all_edges %>%
    mutate(pair = paste(pmin(from, to), pmax(from, to), sep = "|")) %>%
    group_by(pair) %>%
    summarise(count = n(), .groups = 'drop')
  n_boot <- length(cor_list)
  keep_pairs <- edge_counts$pair[edge_counts$count >= n_boot * min_boot]
  bootstrap_net <- all_edges %>%
    mutate(pair = paste(pmin(from, to), pmax(from, to), sep = "|")) %>%
    filter(pair %in% keep_pairs) %>%
    group_by(from, to) %>%
    summarise(rho = mean(rho), fdr = mean(fdr), .groups = 'drop')
  bootstrap_net
}

# Function to calculate node degree and hubs
get_hubs <- function(network, asv_names, n_boot, top_prop = 0.3) {
  # Build igraph object
  g <- graph_from_data_frame(network, directed = FALSE, vertices = asv_names)
  deg <- degree(g, mode = 'all')
  deg_df <- data.frame(node = names(deg), degree = deg)
  deg_df <- deg_df[order(-deg_df$degree), ]
  n_hub <- ceiling(length(deg_df$degree) * top_prop)
  hubs <- deg_df$node[1:n_hub]
  list(degree = deg_df, hubs = hubs)
}

# Main analysis
results <- list()
for (grp in unique(meta$group)) {
  samples_grp <- rownames(meta)[meta$group == grp]
  otu_grp <- otu[, samples_grp, drop=FALSE]
  cor_list <- group_bootstrap_correlation(otu_grp)
  bootstrap_net <- bootstrap_edges(cor_list)
  asv_names <- unique(c(bootstrap_net$from, bootstrap_net$to))
  hub_info <- get_hubs(bootstrap_net, asv_names, length(cor_list))
  results[[grp]] <- list(
    bootstrap_network = bootstrap_net,
    degree = hub_info$degree,
    hubs = hub_info$hubs
  )
}

saveRDS(results, output_file)
cat('Bootstrap co-occurrence network and hubs saved to', output_file, '\n') 