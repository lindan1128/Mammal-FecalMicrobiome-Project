# Cross-disease prediction

Two scripts test whether intestinal disease microbiome signatures generalize beyond their source disease. `cross_disease_prediction.py` asks whether the IBD/diarrhea signature (from LODO meta-analysis) can classify other diseases within the same host, running per-study binary classification in both signature and full-genus modes. `multi_disease_prediction.py` pools canine validation studies and trains a multi-class classifier to distinguish five diseases from Normal simultaneously.

---

## cross_disease_prediction.py

Runs per-study binary classification (sick vs. normal) for 6 bovine and 7 canine validation studies. CZM+CLR transform is applied before splitting; covariate residualization (log sequencing depth + study demographics) and StandardScaler are fit inside each fold. No ComBat: each study is a single batch. Both signature and full genus modes are evaluated.

| Function | Description |
|---|---|
| `load_study_data` | reads `otu_genus_update.csv` and per-study metadata; detects sample ID column by overlap with OTU columns; applies taxonomy aliases for WGS studies |
| `build_feature_matrix` | builds X (samples × genera), y (0/1), and log sequencing depth; missing genera filled with zeros |
| `build_covariate_matrix` | builds numeric covariate matrix from metadata; numeric columns kept as-is, categorical one-hot encoded; returns None if no usable covariates |
| `residualize` | removes covariate effects via OLS residuals; fits on training fold, applies to test |
| `czm_clr` | CZM zero replacement (0.65 × column-min per genus) + CLR per sample; matches lodo_nested scripts |
| `filter_genera_full` | keeps genera with prevalence ≥ 10% and mean relative abundance ≥ 0.01% |
| `safe_auc` / `safe_aupr` | AUROC and AUPR with NaN guard for single-class folds |
| `compute_metrics` | returns AUROC, AUPR, accuracy, sensitivity, specificity, F1 |
| `oversample_balance` | random oversample minority class to match majority class size (for NN) |
| `make_rf` / `make_nn` | base RF and MLP estimators for GridSearchCV |
| `run_cv` | 5-fold stratified CV × 5 repeats; inner 3-fold GridSearchCV; covariate residualization → StandardScaler inside each fold |
| `run_host` | iterates over all validation studies for one host; runs signature and full modes; prints per-study summaries |
| `main` | runs bovine then canine; saves results CSV; prints summary table |

### Key parameters

| Parameter | Bovine | Canine |
|---|---|---|
| Validation studies | 6 (PRJEB25138, PRJNA956982, PRJNA1184454, PRJNA1095683, PRJNA1173238, hy) | 7 (PRJNA474706, PRJNA1012802, PRJNA1167336, PRJNA746550, PRJNA612483, PRJNA837097, PRJNA1150743) |
| Signature genera | 12 (5 enriched, 7 depleted) | 8 (3 enriched, 5 depleted) |
| Outer CV | 5-fold stratified × 5 repeats | 5-fold stratified × 5 repeats |
| Inner CV | 3-fold GridSearchCV | 3-fold GridSearchCV |
| Modes | signature, full | signature, full |
| Covariates | log seq depth + per-study demographics | log seq depth + per-study demographics |

---

## multi_disease_prediction.py

Canine only. Pools 5 validation studies (Parvovirus, CE, Epilepsy, Leishmaniasis, Hyperadrenocorticism) into a single dataset and trains a multi-class classifier predicting disease label or Normal (6 classes total). ComBat batch correction is applied within each CV fold (fit on training split only). Study-level weighting equalizes study contributions for RF; study-aware oversampling does the same for NN.

| Function | Description |
|---|---|
| `load_study_data` | same as cross_disease_prediction.py — reads OTU + metadata, detects sample ID column |
| `czm_clr` | CZM + CLR per sample; applied to full pooled matrix before CV splitting |
| `filter_genera_full` | per-study prevalence + abundance filter |
| `get_common_filtered_genera` | applies per-study filter then intersects across all studies to get full-mode genus list |
| `build_covariate_matrix` | numeric covariate matrix per study; studies without covariates get zero-padded columns |
| `pool_studies` | loads and stacks all canine studies; returns X, disease labels, study batch IDs, log seq depths, padded covariate matrix |
| `apply_combat` | ComBat on combined train+test to avoid pycombat transform() bug; skips if fewer than 2 valid batches |
| `residualize` | OLS covariate residualization inside each fold |
| `compute_study_weights` | per-sample weights so each (class, study) group contributes equally; normalized to balance classes |
| `oversample_study_balanced` | study-aware oversampling for NN: equalizes study sizes within each class, then balances across classes |
| `compute_multiclass_metrics` | weighted AUROC (OVR), weighted AUPR, accuracy, per-class AUROC and AUPR dicts |
| `make_rf` / `make_nn` | base estimators for GridSearchCV |
| `run_multiclass_cv` | 5-fold stratified CV × 5 repeats; covariate residualization → ComBat → StandardScaler → nested 3-fold GridSearchCV inside each fold |
| `main` | runs signature and full modes; saves results CSV; prints per-class summary |

### Key parameters

| Parameter | Value |
|---|---|
| Host | Canine only |
| Studies | 5 (PRJNA474706, PRJNA1012802, PRJNA612483, PRJNA837097, PRJNA1150743) |
| Classes | 6 (Parvovirus, CE, Epilepsy, Leishmaniasis, Hyperadrenocorticism, Normal) |
| Signature genera | 8 (canine LODO signature) |
| Outer CV | 5-fold stratified × 5 repeats |
| Inner CV | 3-fold GridSearchCV (scoring=roc_auc_ovr_weighted) |
| Batch correction | ComBat within each fold (train-only fit) |
| RF class handling | study-level sample weights |
| NN class handling | study-aware oversampling |

---

## Outputs

| File | Contents |
|---|---|
| `tables/cross_disease_prediction_results.csv` | per (host, study, disease, category, mode, model, repeat, fold): AUROC, AUPR, accuracy, sensitivity, specificity, F1 |
| `tables/multi_class_disease_results.csv` | per (mode, model, repeat, fold): weighted AUROC, weighted AUPR, accuracy, per-class AUROC and AUPR columns |
