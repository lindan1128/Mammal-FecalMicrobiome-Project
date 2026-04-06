#!/usr/bin/env python3
# Author: Lin
# Date: 2026-03-20
"""
Cross-disease prediction using bovine diarrhea / canine IBD signature genera.

Tests whether intestinal disease signatures (from LODO meta-analysis) generalize
to other diseases within the same host. For each validation study, applies
CZM+CLR transform and performs 5-fold stratified CV x 5 repeats using RF and NN.
Reports AUROC, AUPR, accuracy, recall, and F1.

Both signature and full genus modes are run. No ComBat is applied because each
validation study is a single study (no multi-study batch effects to correct).
Nested CV: inner 3-fold GridSearchCV for hyperparameter tuning.

Usage:
  python3 cross_disease_prediction.py

Output:
  tables/cross_disease_prediction_results.csv
"""

import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             accuracy_score, recall_score, f1_score)
from sklearn.utils import resample
import warnings
warnings.filterwarnings("ignore")

# config
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = "/Users/lindan/Dropbox/PhD/Projects/Animal_Microbiome_cb"

N_REPEATS = 5
N_FOLDS = 5
RANDOM_STATE = 42

# genus name aliases: older taxonomy (WGS studies) -> current signature name
# applied only when old name exists in OTU table and new name does not
GENUS_ALIASES = {
    "Eubacterium":  "[Eubacterium]_nodatum_group",   # bovine WGS
    "Clostridium":  "Clostridium_sensu_stricto_1",   # canine WGS
    "Prevotella":   "Prevotella_9",                  # canine WGS
}

# signature genera
# bovine: 12 stable genera (proportion >= 3/4 folds), canine: 8 (>= 2/3 folds)
BOVINE_SIGNATURE_GENERA = [
    # Enriched in sick (5)
    "[Ruminococcus]_gauvreauii_group", "[Ruminococcus]_gnavus_group",
    "Lachnospiraceae_UCG-010", "Olsenella", "[Eubacterium]_nodatum_group",
    # Depleted in sick (7)
    "Bacteroides", "Barnesiella", "Clostridia_vadinBB60_group",
    "Muribaculaceae", "Oscillospira", "Rikenellaceae_RC9_gut_group", "Ruminococcus",
]

CANINE_SIGNATURE_GENERA = [
    # Enriched in sick (3)
    "Clostridium_sensu_stricto_1", "Streptococcus", "Terrisporobacter",
    # Depleted in sick (5)
    "Allobaculum", "Megamonas", "Prevotella_9", "Bacteroides", "Peptoclostridium",
]

# validation study configs (disease label, intestinal vs non-intestinal, covariates)
BOVINE_VALIDATION = {
    "PRJEB25138":   {"disease": "Salmonella",     "category": "intestinal",
                     "covariates": ["Farm", "Time"]},
    "PRJNA956982":  {"disease": "Diarrhea",       "category": "intestinal",
                     "covariates": ["Age", "Sex"]},
    "PRJNA1184454": {"disease": "Diarrhea",       "category": "intestinal",
                     "covariates": []},
    "PRJNA1095683": {"disease": "Mastitis",       "category": "non-intestinal",
                     "covariates": []},
    "PRJNA1173238": {"disease": "Low BCS",        "category": "non-intestinal",
                     "covariates": []},
    "hy":           {"disease": "Hyperketonemia", "category": "non-intestinal",
                     "covariates": []},
}

CANINE_VALIDATION = {
    "PRJNA474706":  {"disease": "Parvovirus",          "category": "intestinal",
                     "covariates": []},
    "PRJNA1012802": {"disease": "CE",                  "category": "intestinal",
                     "covariates": []},
    "PRJNA1167336": {"disease": "CE",                  "category": "intestinal",
                     "covariates": []},
    "PRJNA746550":  {"disease": "Epilepsy",            "category": "non-intestinal",
                     "covariates": ["age_(years)", "sex", "BREED"]},
    "PRJNA612483":  {"disease": "Epilepsy",            "category": "non-intestinal",
                     "covariates": ["Age", "Sex", "breed"]},
    "PRJNA837097":  {"disease": "Leishmaniasis",       "category": "non-intestinal",
                     "covariates": ["Age", "Sex", "breed"]},
    "PRJNA1150743": {"disease": "Hyperadrenocorticism","category": "non-intestinal",
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
    """load otu + metadata for one study; handle varied sample ID column names"""
    study_dir = os.path.join(base_dir, study_id)
    otu = pd.read_csv(os.path.join(study_dir, "otu_genus_update.csv"), index_col=0)

    # apply taxonomy aliases: rename old name -> new name when old exists and new does not
    for old, new in GENUS_ALIASES.items():
        if old in otu.index and new not in otu.index:
            otu = otu.rename(index={old: new})

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


def build_feature_matrix(otu, meta, signature_genera):
    """build X (samples x genera) and y (0/1); missing genera filled with zeros"""
    samples = meta["Run"].values
    X = pd.DataFrame(0.0, index=samples, columns=signature_genera)
    avail = [g for g in signature_genera if g in otu.index]
    if avail:
        X[avail] = otu.loc[avail, samples].T.values
    y = (meta.set_index("Run").loc[samples, "group"] == "sick").astype(int).values

    # log sequencing depth per sample
    seq_depth = otu[samples].sum(axis=0).values.astype(float)
    seq_depth[seq_depth == 0] = 1
    log_seq_depth = np.log(seq_depth)

    return X.values, y, log_seq_depth


def build_covariate_matrix(meta, covar_cols):
    """build numeric covariate matrix; numeric columns kept as-is, categorical one-hot encoded"""
    if not covar_cols:
        return None

    samples = meta["Run"].values
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


# transforms and metrics
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
    kept = [
        g for g in otu.index
        if (otu.loc[g] > 0).sum() / n_samples >= min_prevalence
        and relabund.loc[g].mean() >= min_mean_relabund
    ]
    return sorted(kept)


def safe_auc(y_true, y_prob):
    """AUROC; returns NaN if only one class present"""
    if len(np.unique(y_true)) < 2:
        return np.nan
    try:
        return roc_auc_score(y_true, y_prob)
    except ValueError:
        return np.nan


def safe_aupr(y_true, y_prob):
    """AUPR (average precision); returns NaN if only one class"""
    if len(np.unique(y_true)) < 2:
        return np.nan
    try:
        return average_precision_score(y_true, y_prob)
    except ValueError:
        return np.nan


def compute_metrics(y_true, y_pred, y_prob):
    """compute AUROC, AUPR, accuracy, sensitivity, specificity, and F1"""
    auroc = safe_auc(y_true, y_prob)
    aupr = safe_aupr(y_true, y_prob)
    acc = accuracy_score(y_true, y_pred)
    rec = recall_score(y_true, y_pred, zero_division=0)          # sensitivity
    spec = recall_score(y_true, y_pred, pos_label=0, zero_division=0)  # specificity
    f1 = f1_score(y_true, y_pred, zero_division=0)
    return auroc, aupr, acc, rec, spec, f1


def oversample_balance(X, y, random_state=42):
    """random oversample minority class to match majority class size"""
    classes, counts = np.unique(y, return_counts=True)
    max_count = counts.max()
    X_resampled = []
    y_resampled = []
    for cls in classes:
        X_cls = X[y == cls]
        y_cls = y[y == cls]
        if len(X_cls) < max_count:
            X_cls, y_cls = resample(X_cls, y_cls, n_samples=max_count,
                                     replace=True, random_state=random_state)
        X_resampled.append(X_cls)
        y_resampled.append(y_cls)
    return np.vstack(X_resampled), np.concatenate(y_resampled)


# CV pipeline
def run_cv(X, y, study_id, disease, category, host, mode,
           log_seq_depth=None, covar_mat=None):
    """5-fold stratified CV x N repeats for RF and NN on a single study

    CZM+CLR applied before CV; covariate residualization and StandardScaler fit
    inside each fold. Nested CV with inner 3-fold GridSearchCV.
    """
    results = []

    class_counts = np.bincount(y)
    n_min = min(class_counts)
    n_splits = min(N_FOLDS, n_min)
    if n_splits < 2:
        print(f"    WARNING: min class count = {n_min}, skipping CV")
        return results

    X_clr = czm_clr(X)

    # combine log_seq_depth + demographic covariates
    covars_all = None
    if log_seq_depth is not None:
        covars_all = log_seq_depth.reshape(-1, 1)
        if covar_mat is not None:
            covars_all = np.column_stack([covars_all, covar_mat])
    elif covar_mat is not None:
        covars_all = covar_mat

    has_covars = covars_all is not None
    if has_covars:
        print(f"    Covariate adjustment: {covars_all.shape[1]} variables")

    models = {
        "rf": (make_rf, RF_PARAM_GRID),
        "nn": (make_nn, NN_PARAM_GRID),
    }

    for model_name, (model_fn, param_grid) in models.items():
        rskf = RepeatedStratifiedKFold(
            n_splits=n_splits, n_repeats=N_REPEATS, random_state=RANDOM_STATE
        )
        repeat_idx = 0
        fold_in_repeat = 0

        for train_idx, test_idx in rskf.split(X_clr, y):
            X_train, X_test = X_clr[train_idx], X_clr[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            if has_covars:
                X_train, X_test = residualize(
                    X_train, X_test,
                    covars_all[train_idx], covars_all[test_idx])

            scaler = StandardScaler()
            X_train_s = scaler.fit_transform(X_train)
            X_test_s = scaler.transform(X_test)

            # inner 3-fold GridSearchCV for hyperparameter tuning
            base_est = model_fn(random_state=RANDOM_STATE + repeat_idx)
            grid = GridSearchCV(
                base_est, param_grid, cv=3, scoring="roc_auc",
                n_jobs=-1, refit=True,
            )
            # NN: oversample minority class (MLPClassifier has no class_weight)
            if model_name == "nn":
                X_fit, y_fit = oversample_balance(X_train_s, y_train,
                                                   random_state=RANDOM_STATE + repeat_idx)
            else:
                X_fit, y_fit = X_train_s, y_train
            grid.fit(X_fit, y_fit)
            clf = grid.best_estimator_
            print(f"    {model_name.upper()} fold {repeat_idx+1}-{fold_in_repeat+1}: "
                  f"best_params={grid.best_params_}")

            y_pred = clf.predict(X_test_s)
            if hasattr(clf, "predict_proba"):
                y_prob = clf.predict_proba(X_test_s)
                if y_prob.shape[1] == 2:
                    y_prob = y_prob[:, 1]
                else:
                    y_prob = y_prob[:, 0]
            else:
                y_prob = clf.decision_function(X_test_s)

            auroc, aupr, acc, rec, spec, f1 = compute_metrics(y_test, y_pred, y_prob)

            results.append({
                "host": host,
                "study": study_id,
                "disease": disease,
                "category": category,
                "mode": mode,
                "n_genera": X.shape[1],
                "model": model_name,
                "repeat": repeat_idx + 1,
                "fold": fold_in_repeat + 1,
                "auroc": auroc,
                "aupr": aupr,
                "accuracy": acc,
                "sensitivity": rec,
                "specificity": spec,
                "f1": f1,
            })

            fold_in_repeat += 1
            if fold_in_repeat >= n_splits:
                fold_in_repeat = 0
                repeat_idx += 1

    return results


# main pipeline
def run_host(host, validation_config, signature_genera):
    """run cross-disease prediction for one host across all validation studies"""
    if host == "bovine":
        base_dir = os.path.join(PROJECT_DIR, "metadata/bovine")
    else:
        base_dir = os.path.join(PROJECT_DIR, "metadata/canine")

    all_results = []

    for study_id, info in validation_config.items():
        disease = info["disease"]
        category = info["category"]

        print(f"\n  --- {study_id} ({disease}, {category}) ---")

        try:
            otu, meta = load_study_data(base_dir, study_id)
        except Exception as e:
            print(f"    ERROR loading data: {e}")
            continue

        n_sick = sum(meta["group"] == "sick")
        n_normal = sum(meta["group"] == "normal")
        print(f"    Samples: sick={n_sick}, normal={n_normal}, total={n_sick+n_normal}")

        covar_cols = info.get("covariates", [])
        covar_mat = build_covariate_matrix(meta, covar_cols)
        if covar_mat is not None:
            print(f"    Demographic covariates: {covar_cols} -> {covar_mat.shape[1]} features")

        modes = {
            "signature": signature_genera,
            "full": filter_genera_full(otu),
        }
        print(f"    Full mode genera (prev>=10%, relabund>=0.01%): {len(modes['full'])}")

        for mode, genera in modes.items():
            if len(genera) == 0:
                print(f"    WARNING: {mode} mode has 0 genera, skipping")
                continue
            X, y, log_seq_depth = build_feature_matrix(otu, meta, genera)
            print(f"    [{mode}] n_genera={len(genera)}")
            mode_results = run_cv(X, y, study_id, disease, category, host, mode,
                                  log_seq_depth=log_seq_depth, covar_mat=covar_mat)
            all_results.extend(mode_results)

            if mode_results:
                df_tmp = pd.DataFrame(mode_results)
                for model_name in ["rf", "nn"]:
                    sub = df_tmp[df_tmp["model"] == model_name]
                    if len(sub) > 0:
                        auc_vals = sub["auroc"].dropna()
                        aupr_vals = sub["aupr"].dropna()
                        print(f"    {model_name.upper()} [{mode}]: "
                              f"AUC={auc_vals.mean():.3f}+/-{auc_vals.std():.3f}, "
                              f"AUPR={aupr_vals.mean():.3f} (n_folds={len(sub)})")

    return all_results


def main():
    print("=" * 70)
    print("Cross-Disease Prediction Pipeline (nested CV)")
    print("  RF: GridSearchCV(3-fold) over", RF_PARAM_GRID)
    print("  NN: GridSearchCV(3-fold) over", NN_PARAM_GRID)
    print(f"  Outer CV: {N_FOLDS}-fold stratified x {N_REPEATS} repeats")
    print("  Transform: CZM + CLR + covariate residualization")
    print("  Covariates: log_seq_depth + per-study demographics")
    print("  Metrics: AUROC, AUPR, Accuracy, Recall, F1")
    print("=" * 70)

    all_results = []

    # bovine
    print(f"\n{'='*50}")
    print(f"Host: BOVINE ({len(BOVINE_SIGNATURE_GENERA)} signature genera)")
    print(f"{'='*50}")
    bovine_results = run_host("bovine", BOVINE_VALIDATION, BOVINE_SIGNATURE_GENERA)
    all_results.extend(bovine_results)

    # canine
    print(f"\n{'='*50}")
    print(f"Host: CANINE ({len(CANINE_SIGNATURE_GENERA)} signature genera)")
    print(f"{'='*50}")
    canine_results = run_host("canine", CANINE_VALIDATION, CANINE_SIGNATURE_GENERA)
    all_results.extend(canine_results)

    # save
    results_df = pd.DataFrame(all_results)
    out_path = os.path.join(SCRIPT_DIR, "tables", "cross_disease_prediction_results.csv")
    results_df.to_csv(out_path, index=False)
    print(f"\nResults saved to: {out_path}")
    print(f"Total rows: {len(results_df)}")

    # summary
    print("\n" + "=" * 70)
    print("Summary (mean AUC / AUPR +/- sd)")
    print("=" * 70)

    for host_label, host_key in [("BOVINE", "bovine"), ("CANINE", "canine")]:
        print(f"\n{host_label}:")
        df_h = results_df[results_df["host"] == host_key]
        if len(df_h) == 0:
            print("  No results")
            continue

        for model_name in ["rf", "nn"]:
            print(f"  {model_name.upper()}:")
            df_m = df_h[df_h["model"] == model_name]
            for _, grp in df_m.groupby(["study", "disease", "category"]):
                study = grp["study"].iloc[0]
                disease = grp["disease"].iloc[0]
                category = grp["category"].iloc[0]
                auc_vals = grp["auroc"].dropna()
                aupr_vals = grp["aupr"].dropna()
                acc_vals = grp["accuracy"]
                f1_vals = grp["f1"]
                sens_vals = grp["sensitivity"]
                spec_vals = grp["specificity"]
                if len(auc_vals) > 0:
                    print(f"    {study:15s} ({disease:25s}, {category:15s}): "
                          f"AUC={auc_vals.mean():.3f}+/-{auc_vals.std():.3f}, "
                          f"AUPR={aupr_vals.mean():.3f}, "
                          f"Sens={sens_vals.mean():.3f}, "
                          f"Spec={spec_vals.mean():.3f}, "
                          f"F1={f1_vals.mean():.3f}")


if __name__ == "__main__":
    main()
