#!/usr/bin/env Rscript
# Author: Lin
# Date: 2026-03-26
#
# called from lodo_nested_bovine/canine.py for each LODO fold
# does: CZM zero replacement -> CLR -> covariate-adjusted Hedges' g -> DL meta-analysis
#
# usage:
#   Rscript lodo_meta_step1.R <host> <train_studies> <common_genera_file> <output_file> \
#     [p_thresh] [dir_thresh] [min_k]
#
# output: genus, combined_g, p_value, q_value, k, dir_consistency

suppressPackageStartupMessages({
  library(zCompositions)
  library(metafor)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript lodo_meta_step1.R <host> <train_studies> <common_genera_file> <output_file> [p_thresh] [dir_thresh] [min_k]")
}

host              <- args[1]
train_study_ids   <- strsplit(args[2], ",")[[1]]
common_genera_file <- args[3]
output_file       <- args[4]
p_thresh   <- if (length(args) >= 5) as.numeric(args[5]) else 0.1
dir_thresh <- if (length(args) >= 6) as.numeric(args[6]) else 0.5
min_k      <- if (length(args) >= 7) as.integer(args[7]) else 2L

common_genera <- readLines(common_genera_file)
common_genera <- common_genera[nchar(trimws(common_genera)) > 0]

base_dir <- file.path("/Users/lindan/Dropbox/PhD/Projects/Animal_Microbiome_cb/metadata", host)

# covariates per study - only include what's actually available
if (host == "bovine") {
  study_covariates <- list(
    PRJNA716761  = c("Age"),
    PRJNA1263132 = c("Age"),                      # sex all female, breed all dairy calf
    PRJEB25741   = c("Age", "Sex"),               # breed all beef calf
    PRJNA918869  = NULL,                          # no age, sex all female
    PRJNA956982  = c("Age", "Sex")                # breed all beef calf
  )
} else if (host == "canine") {
  study_covariates <- list(
    PRJEB13362  = c("age", "sex", "weight_kg", "breed_grouping"),
    PRJNA905458 = c("Age", "Sex", "Body_weight_(kg)", "Breed"),
    PRJNA861113 = NULL                            # no age/sex, breed all Yorkshire Terrier
  )
} else {
  study_covariates <- list()
}

load_study <- function(study_id) {
  study_dir <- file.path(base_dir, study_id)
  otu  <- read.csv(file.path(study_dir, "otu_genus_update.csv"),
                   row.names = 1, check.names = FALSE)
  meta <- read.csv(file.path(study_dir, paste0(study_id, ".csv")),
                   row.names = 1, check.names = FALSE)
  meta <- meta[meta$group %in% c("normal", "sick"), , drop = FALSE]
  shared <- intersect(colnames(otu), rownames(meta))
  otu  <- otu[, shared, drop = FALSE]
  meta <- meta[shared, , drop = FALSE]
  list(otu = otu, meta = meta, study_id = study_id)
}

compute_clr <- function(study_data, common_genera) {
  otu <- study_data$otu
  mat <- matrix(0, nrow = length(common_genera), ncol = ncol(otu),
                dimnames = list(common_genera, colnames(otu)))
  shared_g <- intersect(common_genera, rownames(otu))
  mat[shared_g, ] <- as.matrix(otu[shared_g, ])

  # filter: >=10% prevalence and >0.01% mean relative abundance
  n_samples <- ncol(mat)
  prevalence <- rowSums(mat > 0) / n_samples
  mat_rel <- sweep(mat, 2, colSums(mat), "/")
  mat_rel[is.nan(mat_rel)] <- 0
  mean_ra <- rowMeans(mat_rel)
  keep <- prevalence >= 0.1 & mean_ra > 0.0001
  mat_detected <- mat[keep, , drop = FALSE]
  detected_genera <- rownames(mat_detected)

  # drop all-zero samples
  sample_ok <- colSums(mat_detected > 0) > 0
  mat_detected <- mat_detected[, sample_ok, drop = FALSE]

  X <- t(mat_detected)  # samples x genera

  # CZM zero replacement; fall back to 0.65 * col-min if cmultRepl errors
  X_rep <- tryCatch({
    suppressWarnings(
      cmultRepl(X, label = 0, method = "CZM", output = "p-counts",
                z.warning = 1.0, z.delete = FALSE)
    )
  }, error = function(e) {
    cat(sprintf("  [CZM fallback] cmultRepl failed for %s: %s\n",
                study_data$study_id, conditionMessage(e)))
    X_fb <- X
    for (j in seq_len(ncol(X_fb))) {
      nz <- X_fb[X_fb[, j] > 0, j]
      if (length(nz) > 0) X_fb[X_fb[, j] == 0, j] <- 0.65 * min(nz)
    }
    X_fb
  })

  X_rep <- as.matrix(X_rep)
  storage.mode(X_rep) <- "double"

  # CLR transform (per sample)
  X_clr <- matrix(NA, nrow = nrow(X_rep), ncol = ncol(X_rep),
                  dimnames = list(rownames(X_rep), colnames(X_rep)))
  for (i in seq_len(nrow(X_clr))) {
    gm <- exp(mean(log(X_rep[i, ])))
    X_clr[i, ] <- log(X_rep[i, ] / gm)
  }

  # fill back into full genus x sample matrix (NA for non-detected genera)
  sample_names <- colnames(otu)
  full_clr <- matrix(NA, nrow = length(common_genera), ncol = length(sample_names),
                     dimnames = list(common_genera, sample_names))
  surviving <- intersect(rownames(X_clr), sample_names)
  full_clr[detected_genera, surviving] <- t(X_clr[surviving, , drop = FALSE])
  full_clr
}

compute_effects_adjusted <- function(clr_mat, meta, study_id, covar_cols, otu_raw) {
  case_ids <- rownames(meta)[meta$group == "sick"]
  ctrl_ids <- rownames(meta)[meta$group == "normal"]
  n1 <- length(case_ids)
  n0 <- length(ctrl_ids)
  genera <- rownames(clr_mat)

  res <- data.frame(genus = genera, study = study_id,
                    n_case = n1, n_ctrl = n0,
                    stringsAsFactors = FALSE)

  sample_ids <- c(case_ids, ctrl_ids)
  group_vec <- c(rep(1, n1), rep(0, n0))
  cov_df <- data.frame(group = group_vec, row.names = sample_ids)

  # log seq depth always included
  seq_depth <- colSums(otu_raw[, sample_ids, drop = FALSE], na.rm = TRUE)
  seq_depth[seq_depth == 0] <- 1
  cov_df$log_seq_depth <- log(seq_depth)
  usable_covars <- "log_seq_depth"

  # add demographic covariates where available and variable
  if (!is.null(covar_cols) && length(covar_cols) > 0) {
    for (cv in covar_cols) {
      if (cv %in% names(meta)) {
        vals <- meta[[cv]][match(sample_ids, rownames(meta))]
        vals_num <- suppressWarnings(as.numeric(vals))
        if (sum(!is.na(vals_num)) > length(sample_ids) * 0.5 &&
            length(unique(vals_num[!is.na(vals_num)])) > 1) {
          cov_df[[make.names(cv)]] <- vals_num
          usable_covars <- c(usable_covars, make.names(cv))
        } else if (is.character(vals) || is.factor(vals)) {
          vals <- as.character(vals)
          vals[vals == "" | is.na(vals)] <- NA
          if (sum(!is.na(vals)) > length(sample_ids) * 0.5 &&
              length(unique(vals[!is.na(vals)])) > 1) {
            cov_df[[make.names(cv)]] <- as.factor(vals)
            usable_covars <- c(usable_covars, make.names(cv))
          }
        }
      }
    }
  }

  for (i in seq_along(genera)) {
    clr_vals <- as.numeric(clr_mat[genera[i], sample_ids])

    if (all(is.na(clr_vals))) {
      res$hedges_g[i] <- res$se_g[i] <- NA
      next
    }

    fit_df <- cov_df
    fit_df$y <- clr_vals
    fit_df <- fit_df[complete.cases(fit_df), ]

    n1_fit <- sum(fit_df$group == 1)
    n0_fit <- sum(fit_df$group == 0)

    if (n1_fit < 3 || n0_fit < 3) {
      res$hedges_g[i] <- res$se_g[i] <- NA
      next
    }

    fml <- as.formula(paste("y ~ group +", paste(usable_covars, collapse = " + ")))
    fit <- tryCatch(lm(fml, data = fit_df), error = function(e) NULL)

    if (!is.null(fit) && "group" %in% rownames(summary(fit)$coefficients)) {
      t_val <- summary(fit)$coefficients["group", "t value"]
      n_total <- n1_fit + n0_fit
      smd <- t_val * sqrt(1/n1_fit + 1/n0_fit)
      J <- 1 - 3 / (4 * (n_total - 2) - 1)
      g_adj <- J * smd
      se_adj <- sqrt(1/n1_fit + 1/n0_fit + g_adj^2 / (2 * n_total))
      res$hedges_g[i] <- g_adj
      res$se_g[i] <- se_adj
    } else {
      # fallback: unadjusted Hedges' g
      x1 <- clr_vals[group_vec == 1]
      x0 <- clr_vals[group_vec == 0]
      md <- mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)
      sd_pool <- sqrt(((n1 - 1) * var(x1, na.rm = TRUE) +
                         (n0 - 1) * var(x0, na.rm = TRUE)) / (n1 + n0 - 2))
      J <- 1 - 3 / (4 * (n1 + n0 - 2) - 1)
      g <- if (sd_pool > 0) J * md / sd_pool else 0
      se_g <- sqrt(1/n1 + 1/n0 + g^2 / (2 * (n1 + n0)))
      res$hedges_g[i] <- g
      res$se_g[i] <- se_g
    }
  }
  res
}

run_meta_and_select <- function(effects_df, p_thresh, dir_thresh, min_k) {
  genera <- unique(effects_df$genus)
  res <- data.frame(genus = genera, stringsAsFactors = FALSE)

  for (i in seq_along(genera)) {
    g <- genera[i]
    df_g <- effects_df[effects_df$genus == g & !is.na(effects_df$se_g) &
                         effects_df$se_g > 0, ]
    k <- nrow(df_g)
    res$k[i] <- k

    if (k < 2) {
      res$combined_g[i] <- res$p_value[i] <- res$dir_consistency[i] <- NA
      next
    }

    tryCatch({
      fit <- rma(yi = df_g$hedges_g, sei = df_g$se_g, method = "DL")
      res$combined_g[i] <- as.numeric(fit$beta)
      res$p_value[i] <- fit$pval
      n_pos <- sum(df_g$hedges_g > 0)
      res$dir_consistency[i] <- max(n_pos, k - n_pos) / k
    }, error = function(e) {
      res$combined_g[i] <<- NA
      res$p_value[i] <<- NA
      res$dir_consistency[i] <<- NA
    })
  }

  # BH correction across all genera
  res$q_value <- p.adjust(res$p_value, method = "BH")

  idx <- !is.na(res$q_value) &
    res$q_value < p_thresh &
    !is.na(res$dir_consistency) &
    res$dir_consistency > dir_thresh &
    res$k >= min_k
  sig <- res[idx, c("genus", "combined_g", "p_value", "q_value", "k", "dir_consistency")]
  sig <- sig[order(sig$q_value), ]

  list(all_meta = res, sig = sig)
}

cat(sprintf("LODO Step 1 [%s]: training on %s\n", host, paste(train_study_ids, collapse = ", ")))
cat(sprintf("  %d common genera | FDR < %g | dir > %g | k >= %d\n",
            length(common_genera), p_thresh, dir_thresh, min_k))

train_data <- setNames(lapply(train_study_ids, load_study), train_study_ids)
clr_data   <- lapply(train_data, compute_clr, common_genera = common_genera)

effects_list <- mapply(
  function(clr_mat, data_obj, study_id) {
    compute_effects_adjusted(clr_mat, data_obj$meta, study_id,
                             study_covariates[[study_id]], data_obj$otu)
  },
  clr_data, train_data, names(train_data),
  SIMPLIFY = FALSE
)
effects_df <- do.call(rbind, effects_list)
rownames(effects_df) <- NULL

result <- run_meta_and_select(effects_df, p_thresh, dir_thresh, min_k)

cat(sprintf("  -> %d signatures (from %d tested)\n",
            nrow(result$sig), nrow(result$all_meta)))
if (nrow(result$sig) > 0) {
  cat("  ", paste(result$sig$genus, collapse = ", "), "\n")
}

write.csv(result$sig, output_file, row.names = FALSE)

full_output <- sub("\\.csv$", "_full.csv", output_file)
write.csv(result$all_meta, full_output, row.names = FALSE)
