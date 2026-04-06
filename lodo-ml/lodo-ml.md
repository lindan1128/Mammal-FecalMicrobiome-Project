# LODO ML pipeline

Five classifiers (LASSO, Ridge, ElasticNet, RF, XGBoost) evaluated under a fully nested LODO scheme. Feature selection is re-run inside each fold using only training studies, so no information from the held-out test study leaks into feature selection.

Two modes per fold: **full** (all common genera) and **signature** (fold-specific genera from nested meta-analysis).

Prerequisite: `lodo_meta_step1.R` must be present in the same directory (called as a subprocess for nested feature selection). See `../meta-analysis/meta-analysis.md`.

---

## lodo_nested_bovine.py / lodo_nested_canine.py

Functions in execution order within `run_lodo()`:

| Function | Description |
|---|---|
| `load_study_data` | reads `otu_genus_update.csv` and per-study metadata; keeps only normal/sick samples and shared sample IDs |
| `select_signature_nested_R` | writes common genera to a temp file, calls `lodo_meta_step1.R` on N−1 training studies (p<0.2, dir_consistency>0.5, k≥2), returns signature genera list and full meta-analysis DataFrame |
| `prepare_Xy` | stacks samples from multiple studies into X/y arrays; applies per-study CZM (0.65 × col-min) + CLR; ensures feature columns are aligned across studies |
| `compute_log2fc_per_study` | CLR-based log2FC (mean_CLR_sick − mean_CLR_normal) / log(2) per genus per study, used for the log2FC heatmap |
| `run_sensitivity_sweep` | 1-D hyperparameter sweep on signature-mode features for each model; each value is evaluated once on the full train→test split without inner CV |
| `fit_model` | inner `StratifiedKFold(n_splits=3)` `GridSearchCV` for hyperparameter tuning; returns test AUC, sensitivity, specificity, best estimator, and best params |

`run_lodo()` ties these together: for each held-out study, calls `select_signature_nested_R` on training folds, prepares X/y in both modes, calls `fit_model` for all five classifiers, extracts feature importance for full mode, and runs the sensitivity sweep on signature features.

### Key parameters

| Parameter | Bovine | Canine |
|---|---|---|
| Studies | PRJNA716761, PRJNA1263132, PRJEB25741, PRJNA918869 | PRJEB13362, PRJNA905458, PRJNA861113 |
| Common genus threshold | ≥3/4 studies | ≥2/3 studies |
| Excluded genera | Dietzia, Corynebacterium, Roseburia, Anaerostipes, Fournierella | — |
| Stability threshold | ≥3/4 folds (ceiling(4 × 3/4)) | ≥2/3 folds (ceiling(3 × 2/3)) |
| Inner CV folds | 3 (reduced if min class < 3) | 3 |

### Model configurations

| Model | Class | Key hyperparameter grid |
|---|---|---|
| LASSO | `LogisticRegression(penalty="l1", solver="saga")` | `C`: [0.001, 0.01, 0.1, 1, 10, 100] |
| Ridge | `LogisticRegression(penalty="l2", solver="saga")` | `C`: [0.001, 0.01, 0.1, 1, 10, 100] |
| ElasticNet | `SGDClassifier(loss="log_loss", penalty="elasticnet")` | `alpha`: [0.0001–1], `l1_ratio`: [0.15, 0.5, 0.85] |
| RF | `RandomForestClassifier(class_weight="balanced")` | `n_estimators`: [50–1000], `max_depth`: [None,10,20], `min_samples_leaf`: [1,3,5] |
| XGBoost | `XGBClassifier(scale_pos_weight=n_neg/n_pos)` | `n_estimators`: [50–500], `max_depth`: [3,6], `learning_rate`: [0.01,0.1] |

### Feature importance

Extracted from `feature_importances_` (RF, XGBoost) or `|coef_|` (LASSO, Ridge, ElasticNet) for full-mode models only. Used by `plot_feature_importance.R` to compare signature vs non-signature genera importance distributions.

---

## Outputs

| File | Contents |
|---|---|
| `bovine/lodo_nested_results.csv` | per (fold, mode, model): AUC, sensitivity, specificity, best_params, n_genera |
| `canine/lodo_nested_results.csv` | same for canine |
| `bovine/lodo_nested_importance.csv` | per (fold, model, genus): raw feature importance, full mode only |
| `canine/lodo_nested_importance.csv` | same for canine |
| `bovine/lodo_hyperparam_sensitivity.csv` | per (fold, model, param_value): AUC from 1-D sweep |
| `canine/lodo_hyperparam_sensitivity.csv` | same for canine |
| `bovine/lodo_nested_fold_signatures.csv` | per-fold genus list; read by `lodo_transfer_learning.py` |
| `canine/lodo_nested_fold_signatures.csv` | same for canine |
| `bovine/lodo_nested_stability.csv` | per genus: times_selected, n_folds, proportion, avg_log2fc |
| `canine/lodo_nested_stability.csv` | same for canine |
| `bovine/lodo_nested_full_meta.csv` | per-fold full meta-analysis results (all genera, not just signatures) |
| `canine/lodo_nested_full_meta.csv` | same for canine |
| `bovine/lodo_nested_log2fc.csv` | log2FC matrix: stable genera × all folds |
| `canine/lodo_nested_log2fc.csv` | same for canine |
