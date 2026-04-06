# Meta-analysis pipeline

Two entry points, both calling `lodo_meta_step1.R` internally:

1. `run_full_meta.R` — full meta-analysis across all studies → final signature genera
2. `lodo_nested_bovine.py` / `lodo_nested_canine.py` — LODO nested CV → stability-selected signatures + ML classification

Final signatures = LODO-stable genera ∩ full meta-analysis genera.

---

## lodo_meta_step1.R

Called as a subprocess by both entry points. Takes a list of training studies and common genera, returns signature genera CSV.

| Function | Description |
|---|---|
| `load_study` | reads OTU table and metadata for one study; keeps only normal/sick samples |
| `compute_clr` | filters genera (≥10% prevalence, >0.01% mean RA), applies CZM zero replacement (`frac=0.65`), then CLR transforms per sample |
| `compute_effects_adjusted` | for each genus, fits OLS regression of CLR values on disease status with log sequencing depth and available demographic covariates; converts t-statistic to bias-corrected Hedges' g |
| `run_meta_and_select` | pools per-study Hedges' g using DL random-effects model (`metafor::rma`), applies BH FDR correction, selects genera meeting FDR < threshold, direction consistency > 0.5, and k ≥ min_k |

---

## run_full_meta.R

No functions. Loops over bovine and canine, finds genera present in ≥3/4 (bovine) or ≥2/3 (canine) studies, writes a temp file, and calls `lodo_meta_step1.R` with all studies as training data.

Outputs: `bovine/full_meta_sig.csv`, `canine/full_meta_sig.csv`

---

## lodo_nested_bovine.py / lodo_nested_canine.py

Functions in execution order within `run_lodo()`:

| Function | Description |
|---|---|
| `load_study_data` | reads OTU table and metadata for one study |
| `select_signature_nested_R` | writes common genera to a temp file, calls `lodo_meta_step1.R` on training studies, returns signature genera and full meta-analysis results |
| `prepare_Xy` | stacks samples from multiple studies into X/y arrays; applies per-study CZM (0.65 × col-min) + CLR for ML features |
| `compute_log2fc_per_study` | computes CLR-based log2FC (sick − normal) per genus per study, used for the heatmap |
| `fit_model` | inner stratified k-fold CV for hyperparameter tuning, returns test AUC |

`run_lodo()` ties these together: for each held-out study, calls `select_signature_nested_R` on training folds, trains classifiers in both signature and full-feature modes, and collects per-fold AUC and log2FC.

A genus is considered stable if selected in ≥3/4 folds (bovine, ≥75%) or ≥2/3 folds (canine, ≥67%).

Outputs: `lodo_nested_results.csv`, `lodo_nested_stability.csv`, `lodo_nested_full_meta.csv`, `lodo_nested_log2fc.csv`

---
