#!/usr/bin/env python3
"""
Cross-disease prediction on other hosts (equine, feline, ovine, swine, human).

Uses bovine (12) + canine (8) signature genera as features.
Same pipeline as cross_disease_prediction.py:
  CZM + CLR -> covariate residualization -> 5-fold stratified CV
  RF (balanced) + NN (oversample) with nested 3-fold GridSearchCV.

Usage:
  python3 cross_disease_other_hosts.py

Output:
  tables/cross_disease_other_hosts_results.csv
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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = "/Users/lindan/Dropbox/PhD/Projects/Animal_Microbiome_cb"

N_REPEATS = 5
N_FOLDS = 5
RANDOM_STATE = 42

# =============================================================================
# Bovine (12) + canine (8) signature genera
# =============================================================================
BOVINE_SIG = [
    # Enriched in sick (2)
    "[Ruminococcus]_gauvreauii_group", "Olsenella",
    # Depleted in sick (10)
    "[Ruminococcus]_torques_group", "Ruminococcus",
    "[Eubacterium]_nodatum_group", "Bacteroides", "Butyricimonas",
    "Barnesiella", "Clostridia_vadinBB60_group",
    "Muribaculaceae", "Oscillospira", "Rikenellaceae_RC9_gut_group",
]
CANINE_SIG = [
    # Enriched in sick (3)
    "Clostridium_sensu_stricto_1", "Streptococcus", "Terrisporobacter",
    # Depleted in sick (5)
    "Bacteroides", "Allobaculum", "Peptoclostridium",
    "Prevotella_9", "Megamonas",
]
SIGNATURE_SETS = {
    "bovine_sig": BOVINE_SIG,
    "canine_sig": CANINE_SIG,
}

# Genus name aliases: older taxonomy -> current signature name.
# Applied when old name is in OTU table and new name is not.
# Needed for WGS studies (PRJNA1218713) and MicrobiomeHD human datasets.
GENUS_ALIASES = {
    "Eubacterium":  "[Eubacterium]_nodatum_group",
    "Clostridium":  "Clostridium_sensu_stricto_1",
    "Prevotella":   "Prevotella_9",
}

# =============================================================================
# Study configurations
# =============================================================================
STUDIES = [
    # --- Equine ---
    {"study": "PRJNA1032075", "host": "equine",
     "disease": "Mixed GI (Horse)",
     "base": os.path.join(PROJECT_DIR, "metadata/equine/PRJNA1032075"),
     "id_col": "Run", "group_col": "group",
     "group_map": None,
     "covariates": ["sex", "Sample_age", "breed_type"]},

    # --- Bovine (Yak) ---
    {"study": "PRJNA1094449", "host": "bovine",
     "disease": "Diarrhea (Yak)",
     "base": os.path.join(PROJECT_DIR, "metadata/bovine/PRJNA1094449"),
     "id_col": "Run", "group_col": "group",
     "group_map": None,
     "covariates": []},
    {"study": "PRJNA761293", "host": "bovine",
     "disease": "Diarrhea (Yak)",
     "base": os.path.join(PROJECT_DIR, "metadata/bovine/PRJNA761293"),
     "id_col": "Run", "group_col": "group",
     "group_map": None,
     "covariates": ["sex"]},

    # --- Feline ---
    {"study": "PRJNA908214", "host": "feline",
     "disease": "IBD (Cat)",
     "base": os.path.join(PROJECT_DIR, "metadata/feline/PRJNA908214"),
     "id_col": "Run", "group_col": "group",
     "group_map": None,
     "covariates": []},
    {"study": "SRP047088", "host": "feline",
     "disease": "CKD (Cat)",
     "base": os.path.join(PROJECT_DIR, "metadata/feline/SRP047088"),
     "id_col": "Run", "group_col": "group",
     "group_map": None,
     "covariates": ["host_sex"]},

    # --- Ovine ---
    {"study": "PRJNA436893", "host": "ovine",
     "disease": "Diarrhea (Sheep)",
     "base": os.path.join(PROJECT_DIR, "metadata/ovine/PRJNA436893"),
     "id_col": "Run", "group_col": "group",
     "group_map": None,
     "covariates": []},

    # --- Swine ---
    {"study": "PRJNA1218713", "host": "swine",
     "disease": "Dysentery (Pig)",
     "base": os.path.join(PROJECT_DIR, "metadata/swine/PRJNA1218713"),
     "id_col": "Run", "group_col": "group",
     "group_map": None,
     "covariates": ["Farm", "Sampling"]},
    {"study": "PRJNA521510", "host": "swine",
     "disease": "Salmonella (Pig)",
     "base": os.path.join(PROJECT_DIR, "metadata/swine/PRJNA521510"),
     "id_col": "Run", "group_col": "group",
     "group_map": None,
     "covariates": ["Pen"],
     "encoding": "latin1"},

    # --- Human (MicrobiomeHD) ---
    {"study": "cdi_schubert", "host": "human",
     "disease": "CDI (Schubert)",
     "base": os.path.join(PROJECT_DIR, "metadata/microbiomehd/cdi_schubert"),
     "id_col": "sample", "group_col": "DiseaseState",
     "group_map": {"CDI": "sick", "H": "normal"},
     "covariates": ["age", "gender"]},
    {"study": "cdi_vincent", "host": "human",
     "disease": "CDI (Vincent)",
     "base": os.path.join(PROJECT_DIR, "metadata/microbiomehd/cdi_vincent"),
     "id_col": "sample", "group_col": "DiseaseState",
     "group_map": {"CDI": "sick", "H": "normal"},
     "covariates": []},
    {"study": "cdi_youngster", "host": "human",
     "disease": "CDI (Youngster)",
     "base": os.path.join(PROJECT_DIR, "metadata/microbiomehd/cdi_youngster"),
     "id_col": "sample", "group_col": "DiseaseState",
     "group_map": {"CDI": "sick", "H": "normal"},
     "covariates": []},
    {"study": "edd_singh", "host": "human",
     "disease": "EDD (Singh)",
     "base": os.path.join(PROJECT_DIR, "metadata/microbiomehd/edd_singh"),
     "id_col": "SampleID", "group_col": "DiseaseState",
     "group_map": {"EDD": "sick", "H": "normal"},
     "covariates": []},
    {"study": "ibd_alm", "host": "human",
     "disease": "IBD (Alm)",
     "base": os.path.join(PROJECT_DIR, "metadata/microbiomehd/ibd_alm"),
     "id_col": "sample", "group_col": "DiseaseState",
     "group_map": {"UC": "sick", "CD": "sick", "nonIBD": "normal"},
     "covariates": ["gender", "age"]},
    {"study": "ibd_engstrand", "host": "human",
     "disease": "IBD (Engstrand)",
     "base": os.path.join(PROJECT_DIR, "metadata/microbiomehd/ibd_engstrand"),
     "id_col": "sample", "group_col": "DiseaseState",
     "group_map": {"UC": "sick", "CD": "sick", "H": "normal"},
     "covariates": ["Gender", "Age"]},
    {"study": "ibd_gevers", "host": "human",
     "disease": "IBD (Gevers)",
     "base": os.path.join(PROJECT_DIR, "metadata/microbiomehd/ibd_gevers"),
     "id_col": "sample", "group_col": "DiseaseState",
     "group_map": {"CD": "sick", "nonIBD": "normal"},
     "covariates": ["race"]},
    {"study": "ibd_huttenhower", "host": "human",
     "disease": "IBD (Huttenhower)",
     "base": os.path.join(PROJECT_DIR, "metadata/microbiomehd/ibd_huttenhower"),
     "id_col": "sample", "group_col": "DiseaseState",
     "group_map": {"UC": "sick", "CD": "sick", "H": "normal"},
     "covariates": ["gender", "age"]},
]

# =============================================================================
# Hyperparameter grids
# =============================================================================
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


# =============================================================================
# Data loading
# =============================================================================
def load_study_data(info):
    enc = info.get("encoding", "utf-8")
    otu = pd.read_csv(os.path.join(info["base"], "otu_genus_update.csv"),
                      index_col=0)

    # Apply taxonomy aliases (old name -> new name, only when old exists and new does not)
    for old, new in GENUS_ALIASES.items():
        if old in otu.index and new not in otu.index:
            otu = otu.rename(index={old: new})

    meta_file = os.path.join(info["base"], f"{info['study']}.csv")
    meta = pd.read_csv(meta_file, encoding=enc)

    id_col = info["id_col"]
    group_col = info["group_col"]

    # Match sample IDs
    otu_samples = set(otu.columns)
    if id_col in meta.columns:
        meta[id_col] = meta[id_col].astype(str)
        shared = list(otu_samples & set(meta[id_col]))
        meta = meta[meta[id_col].isin(shared)].copy()
        meta.index = meta[id_col].values
    else:
        # Try first column or row index
        first_col = meta.columns[0]
        meta[first_col] = meta[first_col].astype(str)
        shared = list(otu_samples & set(meta[first_col]))
        meta = meta[meta[first_col].isin(shared)].copy()
        meta.index = meta[first_col].values

    # Map group labels
    if info["group_map"] is not None:
        raw = meta[group_col].astype(str)
        mapped = raw.map(info["group_map"])
        meta["group"] = mapped
    elif group_col != "group":
        meta["group"] = meta[group_col]

    # Keep only sick/normal
    meta = meta[meta["group"].isin(["sick", "normal"])].copy()
    shared = list(set(otu.columns) & set(meta.index))
    meta = meta.loc[shared]
    otu = otu[shared]

    return otu, meta


def build_feature_matrix(otu, meta, sig_genera):
    samples = list(meta.index)
    X = pd.DataFrame(0.0, index=samples, columns=sig_genera)
    avail = [g for g in sig_genera if g in otu.index]
    if avail:
        X[avail] = otu.loc[avail, samples].T.values
    y = (meta.loc[samples, "group"] == "sick").astype(int).values

    seq_depth = otu[samples].sum(axis=0).values.astype(float)
    seq_depth[seq_depth == 0] = 1
    log_seq_depth = np.log(seq_depth)

    return X.values, y, log_seq_depth


def build_covariate_matrix(meta, covar_cols):
    if not covar_cols:
        return None
    samples = list(meta.index)
    parts = []
    for cv in covar_cols:
        if cv not in meta.columns:
            continue
        vals = meta.loc[samples, cv]
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


# =============================================================================
# CZM + CLR
# =============================================================================
def czm_clr(X):
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


def residualize(X_train, X_test, covars_train, covars_test):
    from numpy.linalg import lstsq
    C_train = np.column_stack([np.ones(len(covars_train)), covars_train])
    C_test = np.column_stack([np.ones(len(covars_test)), covars_test])
    X_train_res, X_test_res = X_train.copy(), X_test.copy()
    for j in range(X_train.shape[1]):
        mask = ~np.isnan(X_train[:, j]) & ~np.any(np.isnan(C_train), axis=1)
        if mask.sum() < C_train.shape[1] + 2:
            continue
        beta, _, _, _ = lstsq(C_train[mask], X_train[mask, j], rcond=None)
        X_train_res[:, j] = X_train[:, j] - C_train @ beta
        X_test_res[:, j] = X_test[:, j] - C_test @ beta
    return X_train_res, X_test_res


# =============================================================================
# Model helpers
# =============================================================================
def make_rf(rs=42):
    return RandomForestClassifier(class_weight="balanced", random_state=rs)

def make_nn(rs=42):
    return MLPClassifier(max_iter=500, early_stopping=False, random_state=rs)

def safe_auc(y_true, y_prob):
    if len(np.unique(y_true)) < 2:
        return np.nan
    try:
        return roc_auc_score(y_true, y_prob)
    except ValueError:
        return np.nan

def safe_aupr(y_true, y_prob):
    if len(np.unique(y_true)) < 2:
        return np.nan
    try:
        return average_precision_score(y_true, y_prob)
    except ValueError:
        return np.nan

def oversample_balance(X, y, rs=42):
    classes, counts = np.unique(y, return_counts=True)
    max_count = counts.max()
    Xs, ys = [], []
    for cls in classes:
        Xc, yc = X[y == cls], y[y == cls]
        if len(Xc) < max_count:
            Xc, yc = resample(Xc, yc, n_samples=max_count, replace=True, random_state=rs)
        Xs.append(Xc); ys.append(yc)
    return np.vstack(Xs), np.concatenate(ys)


# =============================================================================
# CV pipeline
# =============================================================================
def run_cv(X, y, log_seq_depth, covar_mat, info):
    results = []
    class_counts = np.bincount(y)
    n_min = min(class_counts)
    n_splits = min(N_FOLDS, n_min)
    if n_splits < 2:
        print(f"    WARNING: min class={n_min}, skipping")
        return results

    X_clr = czm_clr(X)

    covars_all = log_seq_depth.reshape(-1, 1)
    if covar_mat is not None:
        covars_all = np.column_stack([covars_all, covar_mat])
    has_covars = True
    print(f"    Covariate dims: {covars_all.shape[1]}")

    models = {"rf": (make_rf, RF_PARAM_GRID), "nn": (make_nn, NN_PARAM_GRID)}

    for model_name, (model_fn, param_grid) in models.items():
        rskf = RepeatedStratifiedKFold(
            n_splits=n_splits, n_repeats=N_REPEATS, random_state=RANDOM_STATE)
        rep, fold = 0, 0

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

            base_est = model_fn(rs=RANDOM_STATE + rep)
            grid = GridSearchCV(base_est, param_grid, cv=3,
                                scoring="roc_auc", n_jobs=-1, refit=True)
            if model_name == "nn":
                X_fit, y_fit = oversample_balance(X_train_s, y_train,
                                                  rs=RANDOM_STATE + rep)
            else:
                X_fit, y_fit = X_train_s, y_train
            grid.fit(X_fit, y_fit)
            clf = grid.best_estimator_

            y_pred = clf.predict(X_test_s)
            if hasattr(clf, "predict_proba"):
                y_prob = clf.predict_proba(X_test_s)
                y_prob = y_prob[:, 1] if y_prob.shape[1] == 2 else y_prob[:, 0]
            else:
                y_prob = clf.decision_function(X_test_s)

            results.append({
                "host": info["host"], "study": info["study"],
                "disease": info["disease"],
                "model": model_name,
                "repeat": rep + 1, "fold": fold + 1,
                "auroc": safe_auc(y_test, y_prob),
                "aupr": safe_aupr(y_test, y_prob),
                "accuracy": accuracy_score(y_test, y_pred),
                "sensitivity": recall_score(y_test, y_pred, zero_division=0),
                "specificity": recall_score(y_test, y_pred, pos_label=0, zero_division=0),
                "f1": f1_score(y_test, y_pred, zero_division=0),
            })

            fold += 1
            if fold >= n_splits:
                fold = 0; rep += 1

    return results


# =============================================================================
# Main
# =============================================================================
def main():
    print("=" * 70)
    print("Cross-Disease Prediction — Other Hosts + Human")
    print(f"  Signature sets: bovine_sig ({len(BOVINE_SIG)}), canine_sig ({len(CANINE_SIG)})")
    print(f"  CV: {N_FOLDS}-fold x {N_REPEATS} repeat")
    print(f"  Models: RF + NN (nested 3-fold GridSearchCV)")
    print(f"  Transform: CZM + CLR + covariate residualization")
    print("=" * 70)

    all_results = []

    for info in STUDIES:
        print(f"\n{'='*60}")
        print(f"  {info['study']} ({info['disease']}, {info['host']})")
        print(f"{'='*60}")

        try:
            otu, meta = load_study_data(info)
        except Exception as e:
            print(f"  ERROR loading: {e}")
            continue

        n_sick = (meta["group"] == "sick").sum()
        n_normal = (meta["group"] == "normal").sum()
        print(f"    Samples: sick={n_sick}, normal={n_normal}")

        if n_sick == 0 or n_normal == 0:
            print(f"    SKIP: missing sick or normal samples")
            continue

        covar_mat = build_covariate_matrix(meta, info.get("covariates", []))
        if covar_mat is not None:
            print(f"    Demographic covariates: {info['covariates']} -> {covar_mat.shape[1]} features")

        for sig_name, sig_genera in SIGNATURE_SETS.items():
            print(f"\n  [{sig_name}] ({len(sig_genera)} genera)")
            X, y, log_seq_depth = build_feature_matrix(otu, meta, sig_genera)
            n_avail = sum(1 for g in sig_genera if g in otu.index)
            print(f"    Available: {n_avail}/{len(sig_genera)}")

            study_results = run_cv(X, y, log_seq_depth, covar_mat, info)
            # Tag with signature set
            for r in study_results:
                r["sig_set"] = sig_name
            all_results.extend(study_results)

            if study_results:
                df_tmp = pd.DataFrame(study_results)
                for mn in ["rf", "nn"]:
                    sub = df_tmp[df_tmp["model"] == mn]
                    if len(sub) > 0:
                        auc = sub["auroc"].dropna()
                        aupr = sub["aupr"].dropna()
                        print(f"    {mn.upper()}: AUC={auc.mean():.3f}±{auc.std():.3f}, "
                              f"AUPR={aupr.mean():.3f}, Acc={sub['accuracy'].mean():.3f}")

    # Save
    df = pd.DataFrame(all_results)
    out_path = os.path.join(SCRIPT_DIR, "tables", "cross_disease_other_hosts_results.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")
    print(f"Total rows: {len(df)}")

    # Summary
    print("\n" + "=" * 70)
    print("Summary (RF)")
    print("=" * 70)
    for sig_name in SIGNATURE_SETS:
        print(f"\n  --- {sig_name} ---")
        df_s = df[(df["model"] == "rf") & (df["sig_set"] == sig_name)]
        for _, grp in df_s.groupby(["host", "study", "disease"]):
            r = grp.iloc[0]
            auc = grp["auroc"].dropna()
            aupr = grp["aupr"].dropna()
            print(f"  {r['host']:8s} {r['study']:15s} {r['disease']:25s} "
                  f"AUC={auc.mean():.3f}±{auc.std():.3f}  AUPR={aupr.mean():.3f}")


if __name__ == "__main__":
    main()
