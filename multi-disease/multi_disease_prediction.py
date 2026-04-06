#!/usr/bin/env python3
# Author: Lin
# Date: 2026-03-20
"""
Canine multi-class disease prediction using IBD signature genera.

Pools 5 canine validation studies and trains a multi-class classifier to predict
which disease a sample has. Each sick sample is labeled by disease; all normal
samples are labeled "Normal" (5 diseases + Normal = 6 classes).

Pipeline:
  1. Pool all canine studies -> feature matrix using signature genera
  2. CZM + CLR transform (sample-wise, before splitting)
  3. 5-fold stratified CV x 5 repeats
  4. Within each fold:
     a. ComBat on training fold only (study as batch), transform test fold
     b. StandardScaler: fit on corrected train, transform test
     c. Inner 3-fold GridSearchCV to tune RF or NN hyperparameters
     d. RF: study-level sample weights; NN: study-aware oversampling
     e. Evaluate with weighted AUROC (OVR), weighted AUPR, accuracy, per-class metrics

Usage:
  python3 multi_disease_prediction.py

Output:
  tables/multi_class_disease_results.csv
"""

import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score
from sklearn.utils import resample
from pycombat import Combat
import warnings
warnings.filterwarnings("ignore")

# config
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = "/Users/lindan/Dropbox/PhD/Projects/Animal_Microbiome_cb"

N_REPEATS = 5
N_FOLDS = 5
RANDOM_STATE = 42

# canine signature genera from LODO nested meta-analysis (>= 2/3 folds)
CANINE_SIGNATURE_GENERA = [
    # Enriched in sick (3)
    "Clostridium_sensu_stricto_1", "Streptococcus", "Terrisporobacter",
    # Depleted in sick (5)
    "Allobaculum", "Megamonas", "Prevotella_9", "Bacteroides", "Peptoclostridium",
]

# canine validation studies (IBD meta-analysis studies and WGS studies excluded)
# PRJNA1167336 excluded (WGS shotgun), PRJNA746550 excluded (Epilepsy duplicate)
CANINE_STUDIES = {
    "PRJNA474706":  {"disease": "Parvovirus",          "covariates": []},
    "PRJNA1012802": {"disease": "CE",                  "covariates": []},
    "PRJNA612483":  {"disease": "Epilepsy",
                     "covariates": ["Age", "Sex", "breed"]},
    "PRJNA837097":  {"disease": "Leishmaniasis",
                     "covariates": ["Age", "Sex", "breed"]},
    "PRJNA1150743": {"disease": "Hyperadrenocorticism",
                     "covariates": ["Host_Age", "BREED"]},
}

# hyperparameter grids for nested CV
RF_PARAM_GRID = {
    "n_estimators": [100, 200, 500],
    "max_depth": [None, 10, 20],
    "min_samples_leaf": [1, 3, 5],
}

NN_PARAM_GRID = {
    "hidden_layer_sizes": [(32,), (64, 32), (64, 32, 16)],
    "alpha": [0.0001, 0.001, 0.01],
    "learning_rate_init": [0.001, 0.01],
}


# data loading
def load_study_data(base_dir, study_id):
    """load otu + metadata; handle varied sample ID column names"""
    study_dir = os.path.join(base_dir, study_id)
    otu = pd.read_csv(os.path.join(study_dir, "otu_genus_update.csv"), index_col=0)

    meta = pd.read_csv(os.path.join(study_dir, f"{study_id}.csv"))

    # find which column holds sample IDs matching OTU columns
    otu_samples = set(otu.columns)
    id_col = None

    for candidate in ["Run", "run_accession", "Sample", "Unnamed: 0"]:
        if candidate in meta.columns:
            overlap = len(otu_samples & set(meta[candidate].dropna().astype(str)))
            if overlap > 0:
                id_col = candidate
                break

    # fallback: try the first column (handles BOM-prefixed unnamed columns)
    if id_col is None:
        first_col = meta.columns[0]
        overlap = len(otu_samples & set(meta[first_col].dropna().astype(str)))
        if overlap > 0:
            id_col = first_col

    if id_col is None:
        raise ValueError(
            f"Cannot find sample ID column in {study_id} metadata. "
            f"OTU samples: {list(otu.columns)[:5]}, "
            f"Meta columns: {list(meta.columns)[:10]}"
        )

    if id_col != "Run":
        meta = meta.rename(columns={id_col: "Run"})

    meta = meta[meta["group"].isin(["normal", "sick"])].copy()

    shared = list(set(otu.columns) & set(meta["Run"]))
    if len(shared) == 0:
        raise ValueError(f"No shared samples between OTU and metadata for {study_id}")

    otu = otu[shared]
    meta = meta[meta["Run"].isin(shared)].copy()

    return otu, meta


# transforms
def czm_clr(X):
    """CZM (column-wise, 0.65 x genus-min) + CLR per sample; matches lodo_nested scripts"""
    X = X.copy().astype(float)
    for j in range(X.shape[1]):
        nz = X[X[:, j] > 0, j]
        if len(nz) > 0 and np.any(X[:, j] == 0):
            X[X[:, j] == 0, j] = 0.65 * nz.min()
        elif len(nz) == 0:
            X[:, j] = 1e-10
    log_X = np.log(X)
    gm = log_X.mean(axis=1, keepdims=True)
    clr = log_X - gm
    return np.nan_to_num(clr, nan=0.0, posinf=0.0, neginf=0.0)


def filter_genera_full(otu, min_prevalence=0.10, min_mean_relabund=0.0001):
    """keep genera with prevalence >= 10% and mean relative abundance >= 0.01%"""
    total = otu.sum(axis=0)
    total[total == 0] = 1
    relabund = otu.div(total, axis=1)
    n_samples = otu.shape[1]
    return sorted([
        g for g in otu.index
        if (otu.loc[g] > 0).sum() / n_samples >= min_prevalence
        and relabund.loc[g].mean() >= min_mean_relabund
    ])


def get_common_filtered_genera(host, study_disease_map):
    """per-study prevalence+abundance filter, then intersect across all studies"""
    if host == "bovine":
        base_dir = os.path.join(PROJECT_DIR, "metadata/bovine")
    else:
        base_dir = os.path.join(PROJECT_DIR, "metadata/canine")
    per_study = []
    for study_id in study_disease_map:
        try:
            otu, _ = load_study_data(base_dir, study_id)
            per_study.append(set(filter_genera_full(otu)))
        except Exception as e:
            print(f"  WARNING: Failed to load {study_id} for full-mode filter: {e}")
    if not per_study:
        return []
    common = sorted(set.intersection(*per_study))
    print(f"  Full mode common genera (intersection after per-study filter): {len(common)}")
    return common


# model constructors
def make_rf(random_state=42):
    """base RF estimator for GridSearchCV"""
    return RandomForestClassifier(
        class_weight="balanced",
        random_state=random_state,
    )


def make_nn(random_state=42):
    """base NN estimator for GridSearchCV"""
    return MLPClassifier(
        max_iter=500,
        early_stopping=False,
        random_state=random_state,
    )


# study-level weighting and oversampling
def compute_study_weights(y, batch):
    """per-sample weights so each study contributes equally within each class

    For each (class, study) group: weight = 1 / group_size.
    Normalized within each class so all classes have equal total weight.
    """
    y = np.asarray(y)
    batch = np.asarray(batch)
    weights = np.ones(len(y), dtype=float)

    classes = np.unique(y)
    for cls in classes:
        cls_mask = y == cls
        cls_batch = batch[cls_mask]
        unique_studies = np.unique(cls_batch)
        for study in unique_studies:
            study_mask = cls_batch == study
            n_in_group = study_mask.sum()
            full_idx = np.where(cls_mask)[0][study_mask]
            weights[full_idx] = 1.0 / n_in_group

    # normalize: equal total weight per class
    for cls in classes:
        cls_mask = y == cls
        cls_total = weights[cls_mask].sum()
        if cls_total > 0:
            weights[cls_mask] *= len(y) / (len(classes) * cls_total)

    return weights


def oversample_study_balanced(X, y, batch, random_state=42):
    """study-aware oversampling for NN: equalize study sizes within each class, then balance classes"""
    y = np.asarray(y)
    batch = np.asarray(batch)
    classes = np.unique(y)

    # step 1: within each class, oversample each study to max study size in that class
    class_chunks_X = []
    class_chunks_y = []
    for cls in classes:
        cls_mask = y == cls
        X_cls = X[cls_mask]
        y_cls = y[cls_mask]
        b_cls = batch[cls_mask]
        unique_studies = np.unique(b_cls)
        max_per_study = max((b_cls == s).sum() for s in unique_studies)

        for study in unique_studies:
            s_mask = b_cls == study
            X_s = X_cls[s_mask]
            y_s = y_cls[s_mask]
            if len(X_s) < max_per_study:
                X_s, y_s = resample(X_s, y_s, n_samples=max_per_study,
                                     replace=True, random_state=random_state)
            class_chunks_X.append(X_s)
            class_chunks_y.append(y_s)

    # step 2: balance across classes (oversample to largest class)
    class_sizes = {}
    for cls in classes:
        cls_mask = y == cls
        unique_studies = np.unique(batch[cls_mask])
        max_per_study = max((batch[cls_mask] == s).sum() for s in unique_studies)
        class_sizes[cls] = max_per_study * len(unique_studies)

    max_class_size = max(class_sizes.values())

    final_X = []
    final_y = []
    chunk_idx = 0
    for cls in classes:
        cls_mask = y == cls
        unique_studies = np.unique(batch[cls_mask])
        n_studies = len(unique_studies)
        X_cls_all = np.vstack(class_chunks_X[chunk_idx:chunk_idx + n_studies])
        y_cls_all = np.concatenate(class_chunks_y[chunk_idx:chunk_idx + n_studies])
        chunk_idx += n_studies

        if len(X_cls_all) < max_class_size:
            X_cls_all, y_cls_all = resample(X_cls_all, y_cls_all,
                                             n_samples=max_class_size,
                                             replace=True, random_state=random_state)
        final_X.append(X_cls_all)
        final_y.append(y_cls_all)

    return np.vstack(final_X), np.concatenate(final_y)


# metrics
def compute_multiclass_metrics(y_true, y_pred, y_prob, classes):
    """weighted AUROC (OVR), weighted AUPR, accuracy, and per-class dicts"""
    acc = accuracy_score(y_true, y_pred)

    unique_in_fold = np.unique(y_true)
    if len(unique_in_fold) < 2:
        return np.nan, np.nan, acc, {}, {}

    present_classes = np.unique(y_true)
    try:
        if len(present_classes) == 2:
            pos_cls = present_classes[1]
            auroc_weighted = roc_auc_score(
                (y_true == pos_cls).astype(int), y_prob[:, pos_cls]
            )
        else:
            y_prob_present = y_prob[:, present_classes]
            auroc_weighted = roc_auc_score(
                y_true, y_prob_present, multi_class="ovr", average="weighted",
                labels=present_classes
            )
    except ValueError:
        auroc_weighted = np.nan

    per_class_auroc = {}
    aupr_per_class = []
    weights = []
    for cls in classes:
        y_binary = (y_true == cls).astype(int)
        if y_binary.sum() == 0 or y_binary.sum() == len(y_binary):
            continue
        cls_idx = np.where(classes == cls)[0][0]
        try:
            ap = average_precision_score(y_binary, y_prob[:, cls_idx])
            aupr_per_class.append(ap)
            weights.append(y_binary.sum())
        except ValueError:
            pass
        try:
            from sklearn.metrics import roc_auc_score as _rauc
            auc_cls = _rauc(y_binary, y_prob[:, cls_idx])
            per_class_auroc[int(cls)] = auc_cls
        except ValueError:
            pass

    if len(aupr_per_class) > 0:
        aupr_weighted = np.average(aupr_per_class, weights=weights)
    else:
        aupr_weighted = np.nan

    per_class_aupr = {}
    for cls in classes:
        y_binary = (y_true == cls).astype(int)
        if y_binary.sum() == 0 or y_binary.sum() == len(y_binary):
            continue
        cls_idx = np.where(classes == cls)[0][0]
        try:
            per_class_aupr[int(cls)] = average_precision_score(y_binary, y_prob[:, cls_idx])
        except ValueError:
            pass

    return auroc_weighted, aupr_weighted, acc, per_class_auroc, per_class_aupr


# data pooling
def build_covariate_matrix(meta, covar_cols, samples):
    """build numeric covariate matrix; numeric kept as-is, categorical one-hot encoded"""
    if not covar_cols:
        return None

    parts = []
    for cv in covar_cols:
        if cv not in meta.columns:
            continue
        vals = meta.set_index("Run").loc[samples, cv]
        vals_num = pd.to_numeric(vals, errors="coerce")
        if vals_num.notna().sum() > len(samples) * 0.5 and vals_num.nunique() > 1:
            parts.append(vals_num.fillna(vals_num.median()).values.reshape(-1, 1))
        else:
            vals_cat = vals.astype(str).replace({"nan": np.nan, "": np.nan})
            if vals_cat.notna().sum() > len(samples) * 0.5 and vals_cat.nunique() > 1:
                dummies = pd.get_dummies(vals_cat, drop_first=True).values.astype(float)
                parts.append(dummies)

    if not parts:
        return None
    return np.hstack(parts)


def pool_studies(host, study_disease_map, genera):
    """load and pool all studies; returns X, y_labels, batch, log_seq_depths, covar_mat"""
    if host == "bovine":
        base_dir = os.path.join(PROJECT_DIR, "metadata/bovine")
    else:
        base_dir = os.path.join(PROJECT_DIR, "metadata/canine")

    all_X = []
    all_y = []
    all_batch = []
    all_log_seq = []
    all_covar_parts = []
    max_covar_cols = 0

    # first pass: determine max covariate columns
    study_covars = {}
    for study_id, info in study_disease_map.items():
        covar_cols = info.get("covariates", [])
        try:
            otu, meta = load_study_data(base_dir, study_id)
        except Exception as e:
            print(f"  WARNING: Failed to load {study_id}: {e}")
            continue

        samples = meta["Run"].values
        cmat = build_covariate_matrix(meta, covar_cols, samples)
        if cmat is not None and cmat.shape[1] > max_covar_cols:
            max_covar_cols = cmat.shape[1]
        study_covars[study_id] = (otu, meta, cmat)

    for study_id, info in study_disease_map.items():
        disease_label = info["disease"]
        if study_id not in study_covars:
            continue
        otu, meta, cmat = study_covars[study_id]

        samples = meta["Run"].values
        groups = meta.set_index("Run").loc[samples, "group"].values

        X_study = pd.DataFrame(0.0, index=samples, columns=genera)
        avail = [g for g in genera if g in otu.index]
        if avail:
            X_study[avail] = otu.loc[avail, samples].T.values

        seq_depth = otu[samples].sum(axis=0).values.astype(float)
        seq_depth[seq_depth == 0] = 1
        all_log_seq.append(np.log(seq_depth))

        # pad covariate matrix to max_covar_cols
        if max_covar_cols > 0:
            if cmat is None:
                cmat_padded = np.zeros((len(samples), max_covar_cols))
            elif cmat.shape[1] < max_covar_cols:
                cmat_padded = np.zeros((len(samples), max_covar_cols))
                cmat_padded[:, :cmat.shape[1]] = cmat
            else:
                cmat_padded = cmat
            all_covar_parts.append(cmat_padded)

        for i, sample in enumerate(samples):
            if groups[i] == "sick":
                all_y.append(disease_label)
            else:
                all_y.append("Normal")
            all_batch.append(study_id)

        all_X.append(X_study.values)

        n_sick = (groups == "sick").sum()
        n_normal = (groups == "normal").sum()
        print(f"  {study_id:15s}: sick={n_sick:4d} ({disease_label}), normal={n_normal:4d}")

    X = np.vstack(all_X)
    log_seq_depths = np.concatenate(all_log_seq)
    covar_mat = np.vstack(all_covar_parts) if all_covar_parts else None
    return X, all_y, all_batch, log_seq_depths, covar_mat


# ComBat
def apply_combat(X_train, X_test, batch_train, batch_test):
    """ComBat batch correction within a CV fold; fit on train+test combined to avoid transform() bug"""
    batch_all = list(batch_train) + list(batch_test)
    unique_batches = np.unique(batch_all)
    if len(unique_batches) < 2:
        return X_train, X_test

    # ComBat requires >= 2 samples per batch
    batch_counts = pd.Series(batch_all).value_counts()
    if batch_counts.min() < 2:
        valid_batches = set(batch_counts[batch_counts >= 2].index)
        if len(valid_batches) < 2:
            return X_train, X_test
        mask_train = np.array([b in valid_batches for b in batch_train])
        mask_test = np.array([b in valid_batches for b in batch_test])
        if not mask_train.all() or not mask_test.all():
            return X_train, X_test

    n_train = X_train.shape[0]
    X_combined = np.vstack([X_train, X_test])

    try:
        combat = Combat()
        X_corrected = combat.fit_transform(Y=X_combined, b=batch_all)
        return X_corrected[:n_train], X_corrected[n_train:]
    except Exception as e:
        print(f"    ComBat failed: {e}. Using uncorrected data.")
        return X_train, X_test


# CV pipeline
def residualize(X_train, X_test, covars_train, covars_test):
    """remove covariate effects via OLS residuals; fit on train, apply to test"""
    from numpy.linalg import lstsq

    C_train = np.column_stack([np.ones(len(covars_train)), covars_train])
    C_test = np.column_stack([np.ones(len(covars_test)), covars_test])

    X_train_res = X_train.copy()
    X_test_res = X_test.copy()

    for j in range(X_train.shape[1]):
        mask = ~np.isnan(X_train[:, j]) & ~np.any(np.isnan(C_train), axis=1)
        if mask.sum() < C_train.shape[1] + 2:
            continue
        beta, _, _, _ = lstsq(C_train[mask], X_train[mask, j], rcond=None)
        X_train_res[:, j] = X_train[:, j] - C_train @ beta
        X_test_res[:, j] = X_test[:, j] - C_test @ beta

    return X_train_res, X_test_res


def run_multiclass_cv(X, y_labels, batch_labels, host, label_encoder, mode,
                      log_seq_depths=None, covar_mat=None):
    """5-fold stratified CV x N repeats for canine multi-class prediction

    CZM+CLR before splitting; inside each fold: covariate residualization,
    ComBat (train-only fit), StandardScaler, nested 3-fold GridSearchCV.
    RF uses study-level sample weights; NN uses study-aware oversampling.
    """
    y_encoded = label_encoder.transform(y_labels)
    classes = np.arange(len(label_encoder.classes_))
    batch_arr = np.array(batch_labels)

    # CZM + CLR (per-sample transform before splitting — no leakage)
    X_clr = czm_clr(X)

    covars_all = None
    if log_seq_depths is not None:
        covars_all = log_seq_depths.reshape(-1, 1)
        if covar_mat is not None:
            covars_all = np.column_stack([covars_all, covar_mat])
    elif covar_mat is not None:
        covars_all = covar_mat

    has_covars = covars_all is not None
    if has_covars:
        print(f"  Covariate adjustment: {covars_all.shape[1]} variables")

    results = []
    models = {
        "rf": (make_rf, RF_PARAM_GRID),
        "nn": (make_nn, NN_PARAM_GRID),
    }

    class_counts = np.bincount(y_encoded, minlength=len(classes))
    n_min = class_counts[class_counts > 0].min()
    n_splits = min(N_FOLDS, n_min)
    if n_splits < 2:
        print(f"  WARNING: min class count = {n_min}, cannot do stratified CV")
        return results

    print(f"  Using {n_splits}-fold CV x {N_REPEATS} repeats")

    for model_name, (model_fn, param_grid) in models.items():
        print(f"  Training {model_name.upper()}...")
        rskf = RepeatedStratifiedKFold(
            n_splits=n_splits, n_repeats=N_REPEATS, random_state=RANDOM_STATE
        )
        repeat_idx = 0
        fold_in_repeat = 0

        for train_idx, test_idx in rskf.split(X_clr, y_encoded):
            X_train, X_test = X_clr[train_idx], X_clr[test_idx]
            y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]
            batch_train = batch_arr[train_idx].tolist()
            batch_test = batch_arr[test_idx].tolist()

            if has_covars:
                X_train, X_test = residualize(
                    X_train, X_test,
                    covars_all[train_idx], covars_all[test_idx])

            # ComBat: fit on train only to avoid data leakage
            batch_train_unique = np.unique(batch_train)
            if len(batch_train_unique) >= 2:
                try:
                    combat = Combat()
                    X_train_cb = combat.fit_transform(Y=X_train, b=batch_train)
                    test_in_train = all(b in batch_train_unique for b in batch_test)
                    if test_in_train:
                        X_test_cb = combat.transform(Y=X_test, b=batch_test)
                    else:
                        X_test_cb = X_test
                except Exception:
                    X_train_cb, X_test_cb = X_train, X_test
            else:
                X_train_cb, X_test_cb = X_train, X_test

            scaler = StandardScaler()
            X_train_s = scaler.fit_transform(X_train_cb)
            X_test_s = scaler.transform(X_test_cb)

            unique_train = np.unique(y_train)
            if len(unique_train) < 2:
                fold_in_repeat += 1
                if fold_in_repeat >= n_splits:
                    fold_in_repeat = 0
                    repeat_idx += 1
                continue

            # inner 3-fold GridSearchCV
            base_est = model_fn(random_state=RANDOM_STATE + repeat_idx)
            grid = GridSearchCV(
                base_est, param_grid, cv=3,
                scoring="roc_auc_ovr_weighted",
                n_jobs=-1, refit=True,
            )
            batch_train_arr = np.array(batch_train)
            if model_name == "nn":
                # NN: study-aware oversampling (MLPClassifier has no sample_weight)
                X_fit, y_fit = oversample_study_balanced(
                    X_train_s, y_train, batch_train_arr,
                    random_state=RANDOM_STATE + repeat_idx)
                grid.fit(X_fit, y_fit)
            else:
                # RF: study-level sample weights to equalize study contributions
                sw = compute_study_weights(y_train, batch_train_arr)
                grid.fit(X_train_s, y_train, sample_weight=sw)
            clf = grid.best_estimator_
            print(f"    {model_name.upper()} fold {repeat_idx+1}-{fold_in_repeat+1}: "
                  f"best_params={grid.best_params_}")

            y_pred = clf.predict(X_test_s)
            y_prob = clf.predict_proba(X_test_s)

            # align probability columns with full class set
            y_prob_full = np.zeros((len(y_test), len(classes)))
            for ci, c in enumerate(clf.classes_):
                y_prob_full[:, c] = y_prob[:, ci]

            auroc_w, aupr_w, acc, per_class_auc, per_class_ap = compute_multiclass_metrics(
                y_test, y_pred, y_prob_full, classes
            )

            result = {
                "host": host,
                "mode": mode,
                "n_genera": X.shape[1],
                "model": model_name,
                "repeat": repeat_idx + 1,
                "fold": fold_in_repeat + 1,
                "auroc_weighted": auroc_w,
                "aupr_weighted": aupr_w,
                "accuracy": acc,
            }

            # per-class AUROC and AUPR columns
            for cls_idx, cls_name in enumerate(label_encoder.classes_):
                result[f"auroc_{cls_name}"] = per_class_auc.get(cls_idx, np.nan)
                result[f"aupr_{cls_name}"] = per_class_ap.get(cls_idx, np.nan)

            results.append(result)

            fold_in_repeat += 1
            if fold_in_repeat >= n_splits:
                fold_in_repeat = 0
                repeat_idx += 1

    return results


def main():
    print("=" * 70)
    print("Canine Multi-Class Disease Prediction Pipeline")
    print("  RF: GridSearchCV(3-fold) over", RF_PARAM_GRID)
    print("  NN: GridSearchCV(3-fold) over", NN_PARAM_GRID)
    print(f"  Outer CV: {N_FOLDS}-fold stratified x {N_REPEATS} repeats")
    print("  Inner CV: 3-fold GridSearchCV (scoring=roc_auc_ovr_weighted)")
    print("  Transform: CZM + CLR -> covariate residualization -> ComBat (train-only) -> StandardScaler")
    print("  Metrics: weighted AUROC (OVR), weighted AUPR, accuracy, per-class AUROC/AUPR")
    print("=" * 70)

    # canine multi-class CV
    full_genera = get_common_filtered_genera("canine", CANINE_STUDIES)
    mode_genera = {"signature": CANINE_SIGNATURE_GENERA, "full": full_genera}

    # fit label encoder once (labels independent of genera)
    X_tmp, y_tmp, _, _, _ = pool_studies("canine", CANINE_STUDIES, CANINE_SIGNATURE_GENERA)
    le = LabelEncoder()
    le.fit(sorted(set(y_tmp)))
    print(f"  Classes: {list(le.classes_)}")

    all_results = []
    for mode, genera in mode_genera.items():
        if len(genera) == 0:
            print(f"  [{mode}] WARNING: 0 genera, skipping")
            continue
        print(f"\n  [{mode}] n_genera={len(genera)}")
        X, y, batch, lsd, cov = pool_studies("canine", CANINE_STUDIES, genera)
        results = run_multiclass_cv(X, y, batch, "canine", le, mode,
                                    log_seq_depths=lsd, covar_mat=cov)
        all_results.extend(results)
        if results:
            df_tmp = pd.DataFrame(results)
            for mn in ["rf", "nn"]:
                sub = df_tmp[df_tmp["model"] == mn]
                if len(sub) > 0:
                    auroc = sub["auroc_weighted"].dropna()
                    aupr = sub["aupr_weighted"].dropna()
                    print(f"  {mn.upper()} [{mode}]: "
                          f"AUROC={auroc.mean():.3f}+/-{auroc.std(ddof=0):.3f}, "
                          f"AUPR={aupr.mean():.3f}+/-{aupr.std(ddof=0):.3f}")

    # save
    results_df = pd.DataFrame(all_results)
    out_path = os.path.join(SCRIPT_DIR, "tables", "multi_class_disease_results.csv")
    results_df.to_csv(out_path, index=False)
    print(f"\nResults saved to: {out_path}")
    print(f"Total rows: {len(results_df)}")

    # summary
    print("\n" + "=" * 70)
    print("Summary: Canine Multi-Class CV (mean +/- sd)")
    print("=" * 70)
    if len(results_df) > 0:
        for model_name in ["rf", "nn"]:
            df_m = results_df[results_df["model"] == model_name]
            if len(df_m) == 0:
                continue
            auroc = df_m["auroc_weighted"].dropna()
            aupr = df_m["aupr_weighted"].dropna()
            acc = df_m["accuracy"]
            print(f"  {model_name.upper():3s}: "
                  f"AUROC={auroc.mean():.3f}+/-{auroc.std(ddof=0):.3f}, "
                  f"AUPR={aupr.mean():.3f}+/-{aupr.std(ddof=0):.3f}, "
                  f"Acc={acc.mean():.3f}+/-{acc.std(ddof=0):.3f}")

            auroc_cols = [c for c in df_m.columns if c.startswith("auroc_") and c != "auroc_weighted"]
            if auroc_cols:
                print(f"        Per-class AUROC:")
                for col in sorted(auroc_cols):
                    vals = df_m[col].dropna()
                    cls_name = col.replace("auroc_", "")
                    if len(vals) > 0:
                        print(f"          {cls_name:25s}: {vals.mean():.3f}+/-{vals.std(ddof=0):.3f}")


if __name__ == "__main__":
    main()
