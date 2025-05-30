# Required packages
library(vegan)

otu_file <- 'otu_table.tsv'      # OTU table: rows are ASVs, columns are samples
meta_file <- 'metadata.tsv'      # Metadata: rows are samples, columns include 'group'

# Read data
otu <- read.table(otu_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Ensure sample order consistency
common_samples <- intersect(colnames(otu), rownames(meta))
otu <- otu[, common_samples]
meta <- meta[common_samples, , drop = FALSE]

# Transpose OTU table for vegan (samples as rows)
otu_t <- t(otu)

# PERMANOVA using Bray-Curtis distance
adonis_res <- adonis(otu_t ~ group, data = meta, permutations = 999, method = "bray")

# Print result
print(adonis_res) 