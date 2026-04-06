#!/usr/bin/env python3
# Author: Lin
# Date: 2026-03-25
"""
LODO nested CV for bovine diarrhea.

In each fold, one study is held out; meta-analysis runs on training studies only
(via lodo_meta_step1.R: CZM -> CLR -> covariate-adjusted Hedges' g -> DL model).
Signature genera from training are used as features for the held-out study.

outputs go to bovine/:
  lodo_nested_results.csv, lodo_nested_log2fc.csv,
  lodo_nested_stability.csv, lodo_nested_full_meta.csv
"""

import os
import subprocess
import tempfile
import numpy as np
import pandas as pd
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score
from xgboost import XGBClassifier
import warnings
warnings.filterwarnings("ignore")

# config
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = "/Users/lindan/Dropbox/PhD/Projects/Animal_Microbiome_cb/metadata/bovine"
OUT_DIR  = os.path.join(SCRIPT_DIR, "bovine")
R_SCRIPT = os.path.join(SCRIPT_DIR, "lodo_meta_step1.R")

STUDIES = [
    "PRJNA716761", "PRJNA1263132",
    "PRJEB25741", "PRJNA918869",
    # PRJNA956982 excluded (WGS, not 16S)
    # PRJNA1184454 excluded (WGS, not 16S)
]

PSEUDOCOUNT = 1e-6
P_THRESH = 0.1
DIR_THRESH = 0.5
MIN_K = 2           # LODO has 3 training studies per fold, k>=2 is reasonable
MIN_STUDY_PRESENCE = 3  # genus in >= 3/4 studies
MODEL_CONFIGS = {
    "RF": {
        "estimator": lambda: RandomForestClassifier(class_weight="balanced", random_state=42),
        "param_grid": {
            "n_estimators": [100, 200, 500],
            "max_depth": [None, 10, 20],
            "min_samples_leaf": [1, 3, 5],
        },
    },
    "Ridge": {
        "estimator": lambda: LogisticRegression(penalty="l2", solver="saga",
                                                 class_weight="balanced",
                                                 max_iter=5000, random_state=42),
        "param_grid": {"C": [0.001, 0.01, 0.1, 1, 10]},
    },
    "LASSO": {
        "estimator": lambda: LogisticRegression(penalty="l1", solver="saga",
                                                 class_weight="balanced",
                                                 max_iter=5000, random_state=42),
        "param_grid": {"C": [0.001, 0.01, 0.1, 1, 10]},
    },
    "ElasticNet": {
        "estimator": lambda: SGDClassifier(loss="log_loss", penalty="elasticnet",
                                            class_weight="balanced",
                                            max_iter=5000, random_state=42),
        "param_grid": {"alpha": [0.0001, 0.001, 0.01, 0.1],
                        "l1_ratio": [0.15, 0.5, 0.85]},
    },
    "XGBoost": {
        "estimator": lambda: XGBClassifier(eval_metric="logloss", use_label_encoder=False,
                                            random_state=42),
        "param_grid": {
            "n_estimators": [100, 200],
            "max_depth": [3, 6],
            "learning_rate": [0.01, 0.1],
        },
    },
}

def load_study_data(study_id):
    """load otu + metadata, filter to normal/sick only"""
    study_dir = os.path.join(BASE_DIR, study_id)
    otu  = pd.read_csv(os.path.join(study_dir, "otu_genus_update.csv"), index_col=0)
    meta = pd.read_csv(os.path.join(study_dir, f"{study_id}.csv"))
    if "Run" not in meta.columns and "Unnamed: 0" in meta.columns:
        meta = meta.rename(columns={"Unnamed: 0": "Run"})
    meta = meta[meta["group"].isin(["normal", "sick"])].copy()
    shared = list(set(otu.columns) & set(meta["Run"]))
    otu  = otu[shared]
    meta = meta[meta["Run"].isin(shared)].copy()
    return otu, meta


def select_signature_nested_R(train_studies, common_genera):
    """run lodo_meta_step1.R on training studies, return signature genera + meta results"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        genera_file = f.name
        f.write('\n'.join(common_genera) + '\n')

    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
        output_file = f.name

    try:
        cmd = [
            "Rscript", R_SCRIPT,
            "bovine",
            ",".join(train_studies),
            genera_file,
            output_file,
            str(P_THRESH),
            str(DIR_THRESH),
            str(MIN_K)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(f"  R script stderr: {result.stderr}")
            return [], pd.DataFrame(), pd.DataFrame()

        for line in result.stdout.strip().split('\n'):
            print(f"  [R] {line}")

        full_output = output_file.replace(".csv", "_full.csv")
        sig_genera, sig_df, full_meta_df = [], pd.DataFrame(), pd.DataFrame()
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            sig_df = pd.read_csv(output_file)
            if not sig_df.empty and "genus" in sig_df.columns:
                sig_genera = sorted(sig_df["genus"].tolist())
        if os.path.exists(full_output) and os.path.getsize(full_output) > 0:
            full_meta_df = pd.read_csv(full_output)
        return sig_genera, sig_df, full_meta_df
    finally:
        os.unlink(genera_file)
        if os.path.exists(output_file):
            os.unlink(output_file)
        full_output = output_file.replace(".csv", "_full.csv")
        if os.path.exists(full_output):
            os.unlink(full_output)


def prepare_Xy(study_data, study_ids, genera):
    """stack samples from multiple studies; applies simple CZM + CLR per study"""
    Xs, ys = [], []
    for sid in study_ids:
        otu, meta = study_data[sid]
        samples = meta["Run"].values
        # Build count matrix: samples x genera
        X = pd.DataFrame(0.0, index=samples, columns=genera)
        avail = [g for g in genera if g in otu.index]
        if avail:
            X[avail] = otu.loc[avail, samples].T.values
        mat = X.values.copy()
        # CZM zero replacement: 0.65 * column-min per genus
        for j in range(mat.shape[1]):
            nz = mat[mat[:, j] > 0, j]
            if len(nz) > 0 and np.any(mat[:, j] == 0):
                mat[mat[:, j] == 0, j] = 0.65 * nz.min()
            elif len(nz) == 0:
                mat[:, j] = 1e-10  # genus absent in this study
        # CLR per sample
        log_mat = np.log(mat)
        gm = log_mat.mean(axis=1, keepdims=True)
        clr = log_mat - gm
        # Replace any NaN/Inf with 0
        clr = np.nan_to_num(clr, nan=0.0, posinf=0.0, neginf=0.0)
        y = (meta.set_index("Run").loc[samples, "group"] == "sick").astype(int).values
        Xs.append(clr)
        ys.append(y)
    return np.vstack(Xs), np.concatenate(ys)


def compute_log2fc_per_study(study_data, study_id, genera):
    """log2FC per genus for one study using CLR (sick mean - normal mean)"""
    otu, meta = study_data[study_id]
    avail = [g for g in genera if g in otu.index]
    if not avail:
        return {g: np.nan for g in genera}

    # Count matrix: genera x samples
    mat = otu.loc[avail].values.astype(float)
    sample_ids = list(otu.columns)

    sr = [r for r in meta.loc[meta["group"] == "sick", "Run"].values if r in sample_ids]
    nr = [r for r in meta.loc[meta["group"] == "normal", "Run"].values if r in sample_ids]
    if not sr or not nr:
        return {g: np.nan for g in avail}

    all_samples = sr + nr
    col_idx = [sample_ids.index(s) for s in all_samples]
    X = mat[:, col_idx].T  # samples x genera

    # Remove genera that are all-zero
    detected = X.sum(axis=0) > 0
    X_det = X[:, detected].copy()
    avail_det = [avail[i] for i in range(len(avail)) if detected[i]]

    if X_det.shape[1] == 0:
        return {g: np.nan for g in avail}

    # Zero replacement: 0.65 * column-min (matching R CZM fallback)
    for j in range(X_det.shape[1]):
        nz = X_det[X_det[:, j] > 0, j]
        if len(nz) > 0:
            X_det[X_det[:, j] == 0, j] = 0.65 * nz.min()

    # CLR transform per sample
    log_X = np.log(X_det)
    gm = log_X.mean(axis=1, keepdims=True)
    clr = log_X - gm

    # Mean CLR for sick vs normal
    n_sick = len(sr)
    clr_sick = clr[:n_sick, :].mean(axis=0)
    clr_normal = clr[n_sick:, :].mean(axis=0)

    # Convert CLR difference to log2 scale
    log2fc = (clr_sick - clr_normal) / np.log(2)

    return dict(zip(avail_det, log2fc))


def fit_model(X_train, y_train, X_test, y_test, model_name="RF"):
    """inner-CV hyperparameter tuning + test AUC"""
    cfg = MODEL_CONFIGS[model_name]
    n_splits = min(3, min(np.bincount(y_train)))
    if n_splits < 2:
        n_splits = 2
    inner_cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    est = cfg["estimator"]()
    # Set scale_pos_weight for XGBoost to handle class imbalance
    if model_name == "XGBoost":
        n_neg = np.sum(y_train == 0)
        n_pos = np.sum(y_train == 1)
        est.set_params(scale_pos_weight=n_neg / max(n_pos, 1))
    grid = GridSearchCV(est, cfg["param_grid"], cv=inner_cv, scoring="roc_auc",
                        n_jobs=-1, refit=True, error_score=0.5)
    grid.fit(X_train, y_train)
    if hasattr(grid.best_estimator_, "predict_proba"):
        y_prob = grid.predict_proba(X_test)[:, 1]
    else:
        y_prob = grid.decision_function(X_test)
    try:
        auc = roc_auc_score(y_test, y_prob)
    except ValueError:
        auc = 0.5
    return auc, grid.best_estimator_, grid.best_params_


def run_lodo():
    # load all studies
    study_data = {}
    for sid in STUDIES:
        otu, meta = load_study_data(sid)
        study_data[sid] = (otu, meta)
        print(f"  {sid}: sick={sum(meta['group']=='sick')}, "
              f"normal={sum(meta['group']=='normal')}, genera={otu.shape[0]}")

    # common genera across studies
    genus_counts = Counter()
    for sid in STUDIES:
        genus_counts.update(study_data[sid][0].index.tolist())
    common = sorted([g for g, c in genus_counts.items() if c >= MIN_STUDY_PRESENCE])
    print(f"Common genera (present in >= {MIN_STUDY_PRESENCE}/{len(STUDIES)} studies): {len(common)}")

    results = []
    fold_sig = {}
    fold_meta = {}
    fold_full_meta = {}
    importance_rows = []  # feature importance from full mode

    for test_study in STUDIES:
        train_studies = [s for s in STUDIES if s != test_study]
        print(f"\n--- LODO fold: test = {test_study}, train = {train_studies} ---")

        # nested feature selection via R
        sig_genera, sig_meta, full_meta = select_signature_nested_R(train_studies, common)
        fold_meta[test_study] = sig_meta
        fold_full_meta[test_study] = full_meta
        fold_sig[test_study] = sig_genera
        print(f"  Nested meta-analysis -> {len(sig_genera)} signature genera")

        full_genera = list(common)

        for mode, genera in [("signature", sig_genera), ("full", full_genera)]:
            if not genera:
                print(f"  WARNING: {test_study} [{mode}] -- 0 genera, skipping")
                for mn in MODEL_CONFIGS:
                    results.append(dict(test_study=test_study, mode=mode,
                                        model=mn, n_genera=0, auroc=np.nan,
                                        best_params=""))
                continue

            X_tr, y_tr = prepare_Xy(study_data, train_studies, genera)
            X_te, y_te = prepare_Xy(study_data, [test_study], genera)

            for model_name in MODEL_CONFIGS:
                auc, clf, best_params = fit_model(
                    X_tr, y_tr, X_te, y_te, model_name)
                results.append(dict(test_study=test_study, mode=mode,
                                    model=model_name, n_genera=len(genera),
                                    auroc=auc, best_params=str(best_params)))
                print(f"  [{mode:9s}] {model_name:5s}: AUC={auc:.3f}  n={len(genera)}  {best_params}")

                # Extract feature importance for full mode
                if mode == "full":
                    imp = None
                    if hasattr(clf, "feature_importances_"):
                        imp = clf.feature_importances_
                    elif hasattr(clf, "coef_"):
                        imp = np.abs(clf.coef_).flatten()
                    if imp is not None and len(imp) == len(genera):
                        for g, v in zip(genera, imp):
                            importance_rows.append(dict(
                                fold_test=test_study, model=model_name,
                                genus=g, importance=v
                            ))

    # stability: selected in >= min_stable folds
    min_stable = int(np.ceil(len(STUDIES) * 3 / 4))
    genus_count = {}
    for sigs in fold_sig.values():
        for g in sigs:
            genus_count[g] = genus_count.get(g, 0) + 1
    stable_genera = sorted([g for g, c in genus_count.items() if c >= min_stable])
    print(f"\nStable genera (>={min_stable}/{len(STUDIES)} folds): {len(stable_genera)}")

    # log2FC per study for stable genera
    study_log2fc = {}
    for sid in STUDIES:
        study_log2fc[sid] = compute_log2fc_per_study(
            study_data, sid, stable_genera)

    # Save feature importance
    if importance_rows:
        imp_df = pd.DataFrame(importance_rows)
        imp_df.to_csv(os.path.join(OUT_DIR, "lodo_nested_importance.csv"), index=False)
        print(f"Feature importance saved: {len(imp_df)} rows")

    return pd.DataFrame(results), study_log2fc, fold_meta, fold_sig, fold_full_meta



if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"LODO bovine | output: {OUT_DIR}")

    results_df, fold_log2fc, fold_meta, fold_sig, fold_full_meta = run_lodo()

    # save outputs
    results_df.to_csv(os.path.join(OUT_DIR, "lodo_nested_results.csv"), index=False)

    # per-fold full meta results
    full_meta_rows = []
    for test_study, fm in fold_full_meta.items():
        if not fm.empty:
            fm_copy = fm.copy()
            fm_copy["fold_test"] = test_study
            full_meta_rows.append(fm_copy)
    if full_meta_rows:
        pd.concat(full_meta_rows, ignore_index=True).to_csv(
            os.path.join(OUT_DIR, "lodo_nested_full_meta.csv"), index=False)

    # log2FC matrix
    all_g = sorted(set().union(*[d.keys() for d in fold_log2fc.values()]))
    coef_df = pd.DataFrame(index=all_g, columns=STUDIES)
    for sid, d in fold_log2fc.items():
        for g, w in d.items():
            coef_df.loc[g, sid] = w
    coef_df.index.name = "genus"
    coef_df.to_csv(os.path.join(OUT_DIR, "lodo_nested_log2fc.csv"))

    # stability summary
    genus_count = {}
    for sigs in fold_sig.values():
        for g in sigs:
            genus_count[g] = genus_count.get(g, 0) + 1
    avg_log2fc_all = {}
    for g in genus_count:
        vals = [d.get(g) for d in fold_log2fc.values() if g in d]
        avg_log2fc_all[g] = np.mean(vals) if vals else np.nan
    stab = pd.DataFrame([
        {"genus": g, "times_selected": c, "n_folds": len(STUDIES),
         "proportion": c / len(STUDIES),
         "avg_log2fc": avg_log2fc_all.get(g, np.nan)}
        for g, c in sorted(genus_count.items(), key=lambda x: -x[1])
    ])
    stab.to_csv(os.path.join(OUT_DIR, "lodo_nested_stability.csv"), index=False)

    # Save per-fold signatures
    fold_sig_rows = []
    for sid, sigs in fold_sig.items():
        for g in sigs:
            fold_sig_rows.append({"fold_test": sid, "genus": g})
    if fold_sig_rows:
        pd.DataFrame(fold_sig_rows).to_csv(
            os.path.join(OUT_DIR, "lodo_nested_fold_signatures.csv"), index=False)
        print(f"Per-fold signatures saved: {len(fold_sig_rows)} rows")

    print()
    for mn in MODEL_CONFIGS:
        for mode in ["full", "signature"]:
            a = results_df[(results_df["model"] == mn) &
                           (results_df["mode"] == mode)]["auroc"].dropna()
            if len(a) > 0:
                print(f"  {mn:5s} {mode:10s}: mean AUC = {a.mean():.3f} +/- {a.std():.3f}")

    min_stable = int(np.ceil(len(STUDIES) * 3 / 4))
    print(f"\nStable signature (>={min_stable}/{len(STUDIES)} folds):")
    for g in all_g:
        fc = coef_df.loc[g].astype(float).mean()
        d = "sick-enriched" if fc > 0 else "normal-enriched"
        print(f"  {g:45s} {d}  avg log2FC={fc:.3f}")
