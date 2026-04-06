#!/usr/bin/env Rscript
# run full meta-analysis (all studies) to get final signature genera
# outputs go to bovine/ and canine/ subdirectories

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
base_dir   <- "/Users/lindan/Dropbox/PhD/Projects/Animal_Microbiome_cb/metadata"

hosts <- list(
  bovine = list(
    studies = c("PRJNA716761", "PRJNA1263132",
                "PRJEB25741", "PRJNA918869"),
    # PRJNA956982 and PRJNA1184454 excluded (WGS, not 16S)
    min_presence = 3,   # genus must appear in >= 3/4 studies
    min_k = 3
  ),
  canine = list(
    studies = c("PRJEB13362", "PRJNA905458", "PRJNA861113"),
    min_presence = 2,
    min_k = 3           # all 3 studies must contribute
  )
)

p_thresh   <- 0.1  # FDR threshold (BH correction applied in lodo_meta_step1.R)
dir_thresh <- 0.5

for (host in names(hosts)) {
  cfg     <- hosts[[host]]
  studies <- cfg$studies
  cat(sprintf("\n--- %s (%d studies) ---\n", host, length(studies)))

  # collect genera present in each study
  genera_per_study <- list()
  for (sid in studies) {
    otu_file <- file.path(base_dir, host, sid, "otu_genus_update.csv")
    otu_full <- read.csv(otu_file, row.names = 1, check.names = FALSE)
    genera_per_study[[sid]] <- rownames(otu_full)
    cat(sprintf("  %s: %d genera\n", sid, length(genera_per_study[[sid]])))
  }

  # keep genera present in >= min_presence studies
  all_genera <- unique(unlist(genera_per_study))
  presence   <- sapply(all_genera, function(g) {
    sum(sapply(genera_per_study, function(x) g %in% x))
  })
  common <- all_genera[presence >= cfg$min_presence]
  cat(sprintf("  common genera (>= %d/%d studies): %d\n",
              cfg$min_presence, length(studies), length(common)))

  tmp_genera <- tempfile(fileext = ".txt")
  writeLines(common, tmp_genera)

  out_dir <- file.path(script_dir, host)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_file <- file.path(out_dir, "full_meta_sig.csv")

  r_script <- file.path(script_dir, "lodo_meta_step1.R")
  cmd <- sprintf(
    'Rscript "%s" %s %s "%s" "%s" %g %g %d',
    r_script, host, paste(studies, collapse = ","),
    tmp_genera, out_file, p_thresh, dir_thresh, cfg$min_k
  )
  system(cmd)

  if (file.exists(out_file)) {
    sig <- read.csv(out_file)
    cat(sprintf("\n  %s: %d signature genera\n", host, nrow(sig)))
    if (nrow(sig) > 0) {
      enriched <- sig$genus[sig$combined_g > 0]
      depleted <- sig$genus[sig$combined_g < 0]
      cat(sprintf("    enriched (%d): %s\n", length(enriched), paste(enriched, collapse = ", ")))
      cat(sprintf("    depleted (%d): %s\n", length(depleted), paste(depleted, collapse = ", ")))
    }
  }

  unlink(tmp_genera)
}

cat("\ndone.\n")
