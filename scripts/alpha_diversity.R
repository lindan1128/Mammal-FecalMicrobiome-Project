# Required packages
library(phyloseq)
library(dplyr)

otu_file <- 'otu_table.tsv'      # OTU table: rows are ASVs, columns are samples
meta_file <- 'metadata.tsv'      # Metadata: rows are samples, columns include 'group'

# Read data
otu <- read.table(otu_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Create phyloseq object
OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
SAMPLEDATA <- sample_data(meta)
physeq_sub <- phyloseq(OTU, SAMPLEDATA)

# Calculate alpha diversity for each group
shannon_normal <- estimate_richness(subset_samples(physeq_sub, group == "normal"), measures = "Shannon")
shannon_disease <- estimate_richness(subset_samples(physeq_sub, group == "disease"), measures = "Shannon")
chao1_normal <- estimate_richness(subset_samples(physeq_sub, group == "normal"), measures = "Chao1")
chao1_disease <- estimate_richness(subset_samples(physeq_sub, group == "disease"), measures = "Chao1")

# Combine results
shannon_df <- rbind(
  data.frame(group = "normal", value = shannon_normal$Shannon),
  data.frame(group = "disease", value = shannon_disease$Shannon)
)
chao1_df <- rbind(
  data.frame(group = "normal", value = chao1_normal$Chao1),
  data.frame(group = "disease", value = chao1_disease$Chao1)
)

# Wilcoxon rank-sum test
shannon_p <- wilcox.test(value ~ group, data = shannon_df)$p.value
chao1_p <- wilcox.test(value ~ group, data = chao1_df)$p.value

# Summary statistics
summary_stats <- data.frame(
  group = c("normal", "disease"),
  Shannon_mean = c(mean(shannon_normal$Shannon, na.rm = TRUE), mean(shannon_disease$Shannon, na.rm = TRUE)),
  Shannon_sd = c(sd(shannon_normal$Shannon, na.rm = TRUE), sd(shannon_disease$Shannon, na.rm = TRUE)),
  Chao1_mean = c(mean(chao1_normal$Chao1, na.rm = TRUE), mean(chao1_disease$Chao1, na.rm = TRUE)),
  Chao1_sd = c(sd(chao1_normal$Chao1, na.rm = TRUE), sd(chao1_disease$Chao1, na.rm = TRUE))
)

# Return results as a list
result <- list(
  summary_stats = summary_stats,
  shannon_p = shannon_p,
  chao1_p = chao1_p
)
result 