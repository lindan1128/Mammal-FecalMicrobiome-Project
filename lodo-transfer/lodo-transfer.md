# Transfer learning pipeline

Neural network (MLP) transfer learning evaluated under the same LODO folds as the ML pipeline. Three strategies are compared on the same held-out 25% test split within each of 20 paired repeats: LODO cross-validation only, transfer (fine-tuning on target domain), and transfer with ComBat batch correction on the source domain.

Two feature sets per fold: **signature** (per-fold genera from `lodo_nested_fold_signatures.csv`, no leakage) and **full** (all common genera matching `lodo_nested_*.py`).

Prerequisite: `lodo_nested_bovine.py` / `lodo_nested_canine.py` must have been run first to produce `lodo_nested_fold_signatures.csv`.

---

## lodo_transfer_learning.py

### Top-level functions

| Function | Description |
|---|---|
| `load_study_data` | reads OTU table and metadata for one study; keeps only normal/sick samples and shared sample IDs |
| `prepare_Xy_with_batch` | stacks samples from multiple studies; applies per-study CZM (0.65 × col-min) + CLR; also returns a batch label array (study ID per sample) for ComBat |
| `get_fold_signatures` | reads `lodo_nested_fold_signatures.csv` and returns a dict `{test_study: [genus, ...]}` |
| `get_common_genera` | builds the full common-genus list matching `lodo_nested_*.py` (presence threshold + exclusions) |
| `apply_combat` | fits `pycombat.Combat` on source training data using study ID as batch; returns corrected X; falls back to uncorrected if <2 batches or Combat fails |
| `make_nn` | constructs `MLPClassifier` with fixed alpha=0.001, lr=0.001, early stopping (validation_fraction=0.15, n_iter_no_change=20) |
| `safe_metrics` | computes AUC, sensitivity, and specificity from probability scores; returns NaN safely if only one class present |
| `cv_select_architecture` | inner `StratifiedKFold(n_splits=3)` over the 3 candidate architectures on source training data; returns the architecture with highest mean AUC |
| `paired_eval` | runs all three strategies across N_REPEATS=20 `StratifiedShuffleSplit(test_size=0.25)` splits; all three strategies share the same 25% test indices within each repeat |
| `arch_eval` | evaluates all 3 architectures under the transfer (no ComBat) strategy across the same 20 splits; no inner CV selection, used for architecture comparison figure |
| `run_host` | orchestrates all folds and feature sets for one host; calls `cv_select_architecture`, `paired_eval`, and `arch_eval` (signature only); collects results and architecture rows |

### Three strategies in `paired_eval`

| Strategy | Source training | Target use | ComBat |
|---|---|---|---|
| LODO | N−1 studies (all samples) → test 25% directly | none | no |
| Transfer | N−1 studies pretrain → fine-tune on 75% → test 25% | warm_start, max_iter=200 | no |
| Transfer+ComBat | ComBat-corrected N−1 pretrain → fine-tune on 75% → test 25% | warm_start, max_iter=200 | source only |

ComBat is fit and applied to the source (N−1) training matrix only. The 75% fine-tune set and 25% test set are never passed through ComBat.

### Architectures

| Name | hidden_layer_sizes |
|---|---|
| shallow | (32,) |
| medium | (64, 32) |
| deep | (64, 32, 16) |

Architecture for LODO and Transfer strategies is selected via `cv_select_architecture` on raw source data. Transfer+ComBat selects architecture separately on ComBat-corrected source data.

### Key parameters

| Parameter | Value |
|---|---|
| Repeats | 20 (`StratifiedShuffleSplit`, test_size=0.25, random_state=42) |
| Fixed alpha | 0.001 |
| Fixed learning rate | 0.001 |
| Pretrain max_iter | 500 |
| Fine-tune max_iter | 200 (warm_start=True) |
| Inner CV folds (arch selection) | 3 |

---

## Outputs

| File | Contents |
|---|---|
| `tables/transfer_learning_comprehensive_results.csv` | one row per (host, test_study, feature_set, strategy, combat, repeat): AUC, sensitivity, specificity, best_architecture, n_genera |
| `tables/transfer_architecture_results.csv` | one row per (host, test_study, architecture, repeat): AUC, sensitivity, specificity; signature feature set, transfer strategy only |
