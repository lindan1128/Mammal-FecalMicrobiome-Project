### Toward diseases-associated fecal microbiome signatures for domestic mammals

Reproducibility scripts for the analysis pipeline. Each subdirectory contains scripts and a detailed method description for one analytical module.

<img src="https://github.com/lindan1128/Mammal-FecalMicrobiome-Project/blob/main/Fig1.jpg" alt="Workflow diagram">

```
github/
├── meta-analysis/                    # LODO nested meta-analysis to derive signature genera
│   ├── meta-analysis.md              # Method description and parameter reference
│   ├── run_full_meta.R               # Full meta-analysis across all studies
│   ├── lodo_meta_step1.R             # Core meta-analysis subroutine (called per fold)
│   ├── lodo_nested_bovine.py         # LODO nested CV — bovine (4 studies)
│   └── lodo_nested_canine.py         # LODO nested CV — canine (3 studies)
├── lodo-ml/                          # LODO machine learning
│   ├── lodo-ml.md                    # Method description and parameter reference
│   ├── lodo_nested_bovine.py         # LODO ML pipeline — bovine
│   └── lodo_nested_canine.py         # LODO ML pipeline — canine
├── lodo-transfer/                    # LODO transfer learning 
│   ├── lodo-transfer.md              # Method description and parameter reference
│   └── lodo_transfer_learning.py     # LODO CV vs. LODO Transfer vs. LODO Transfer+ComBat
├── multi-disease/                    # Signature generalisation within host species
│   ├── cross-diseases.md             # Method description and parameter reference
│   ├── cross_disease_prediction.py   # Cross-disease binary classification
│   └── multi_disease_prediction.py   # Multi-class disease classification (canine)
├── cross-host/                       # Signature generalisation across host species
│   └── cross_disease_other_hosts.py  # Cross-host prediction (equine, caprine, feline, swine, human)
├── tables/                           # Source data for main figures (one .xlsx per panel)
│   ├── fig2c.xlsx                    # Fig. 2c — phylum-level relative abundance
│   ├── fig2de.xlsx                   # Fig. 2d–e — PERMANOVA R² (Bray-Curtis) and Shannon log2FC
│   ├── fig2fghi.xlsx                 # Fig. 2f–i — MDS and dbRDA coordinates, PERMANOVA tables
│   ├── fig3a.xlsx                    # Fig. 3a — LODO forest plot (Hedges' g) and log2FC heatmap
│   ├── fig3b.xlsx                    # Fig. 3b — LODO AUC per study (post-adjustment)
│   ├── fig3c.xlsx                    # Fig. 3c — Transfer learning AUC comparison (post-adjustment)
│   ├── fig4ab.xlsx                   # Fig. 4a–b — Cross-disease AUC, bovine and canine (post-adjustment)
│   ├── fig4cd.xlsx                   # Fig. 4c–d — Bray-Curtis R² vs. cross-disease AUC
│   ├── fig4ef.xlsx                   # Fig. 4e–f — Signature heatmap log2FC scores, bovine and canine
│   ├── fig4g.xlsx                    # Fig. 4g — Multi-class AUC per disease class (canine)
│   ├── fig5b.xlsx                    # Fig. 5b — Cross-host AUC, bovine and canine signatures (post-adjustment)
│   ├── fig5c.xlsx                    # Fig. 5c — Other-hosts signature heatmap log2FC scores
│   └── fig5d.xlsx                    # Fig. 5d — Signature genera log2(case/control) across hosts
└── README.md
```

---

## 1. meta-analysis/

Derives disease-associated signature genera via Leave-One-Dataset-Out (LODO) nested meta-analysis using a DerSimonian–Laird random-effects model with Hedges' *g* effect sizes.

| Script | Description |
|--------|-------------|
| `lodo_meta_step1.R` | Core meta-analysis subroutine: CZM zero replacement → CLR transformation → covariate-adjusted Hedges' *g* → DL random-effects model → BH FDR correction. Called by both `run_full_meta.R` and the LODO Python scripts. |
| `lodo_nested_bovine.py` | LODO nested CV for bovine (4 studies): in each fold, meta-analysis is re-run on the 3 training studies, and signature genera are used to classify the held-out study. |
| `lodo_nested_canine.py` | Same as above for canine (3 studies). |

**Key outputs:** `lodo_nested_results.csv`, `lodo_nested_stability.csv`, `lodo_nested_log2fc.csv`, `lodo_nested_full_meta.csv`

---

## 2. lodo-ml/

Evaluates multiple classifiers under the LODO nested CV framework. Feature selection is performed inside each fold to prevent data leakage.

| Script | Description |
|--------|-------------|
| `lodo_nested_bovine.py` | LODO ML pipeline for bovine: five classifiers (LASSO, Ridge, ElasticNet, Random Forest, XGBoost) in both signature and full-feature modes; hyperparameter sensitivity sweep; feature importance extraction. |
| `lodo_nested_canine.py` | Same pipeline for canine. |

**Key outputs:** `lodo_nested_results.csv`, `lodo_nested_stability.csv`, `lodo_nested_importance.csv`, `lodo_hyperparam_sensitivity.csv`, `lodo_nested_fold_signatures.csv`

---

## 3. lodo-transfer/

Compares three transfer learning strategies on identical held-out test sets, using per-fold signature genera from the LODO meta-analysis (no leakage).

| Script | Description |
|--------|-------------|
| `lodo_transfer_learning.py` | Evaluates LODO (train on N−1 studies, test directly), Transfer (pretrain on N−1, fine-tune on 75% target), and Transfer+ComBat (ComBat batch correction on source before pretraining) across 20 paired repeats and three MLP architectures. Architecture selected via inner 3-fold CV. Learning rate: 0.0001. |

**Key outputs:** `transfer_learning_comprehensive_results.csv`, `transfer_architecture_results.csv`

---

## 4. multi-disease/

Applies the disease signatures to multi-disease classification scenarios within each host species.

| Script | Description |
|--------|-------------|
| `cross_disease_prediction.py` | Tests whether bovine/canine diarrhea signatures generalise to other diseases within the same host. Binary classification (sick vs. normal) across 6 bovine + 7 canine validation studies; 5-fold stratified CV × 5 repeats; Random Forest and Neural Network. |
| `multi_disease_prediction.py` | Multi-class disease classification in canine: pools 5 validation studies to predict 6 classes (5 diseases + normal). ComBat batch correction applied within each CV fold. Signature and full-feature modes. |

**Key outputs:** `cross_disease_prediction_results.csv`, `multi_class_disease_results.csv`

---

## 5. cross-host/

Tests whether bovine and canine signatures generalise to other host species and to human gut microbiome studies.

| Script | Description |
|--------|-------------|
| `cross_disease_other_hosts.py` | Applies the bovine (12 genera) and canine (8 genera) signatures to equine, feline, ovine, swine, and human datasets. CZM+CLR preprocessing, covariate residualisation, and 5-fold stratified CV × 5 repeats with Random Forest and Neural Network. |

**Key outputs:** `cross_disease_other_hosts_results.csv`

---

## Dependencies

**R:** `meta`, `vegan`, `dplyr`, `ggplot2`, `svglite`, `patchwork`, `sva`

**Python:** `scikit-learn`, `numpy`, `pandas`, `scipy`, `neurocombat_sklearn`
