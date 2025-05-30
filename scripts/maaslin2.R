# Required packages
library(Maaslin2)
library(dplyr)

otu_file <- 'otu_table.tsv'      # OTU table: rows are taxa, columns are samples
meta_file <- 'metadata.tsv'      # Metadata: rows are samples, columns include 'group'

# Read data
otu <- read.table(otu_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Ensure sample order consistency
common_samples <- intersect(colnames(otu), rownames(meta))
otu <- otu[, common_samples]
meta <- meta[common_samples, , drop = FALSE]

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = t(otu),                # MaAsLin2 expects samples as rows
  input_metadata = meta,
  output = NULL,
  fixed_effects = "group",
  reference = list(group = "normal"),
  normalization = "NONE",             # Set as needed, e.g., "TSS" for relative abundance
  transform = "NONE"                  # Set as needed, e.g., "LOG"
)

# Print significant results (FDR < 0.1)
fit_data$results %>% filter(qval < 0.1) %>% print() 