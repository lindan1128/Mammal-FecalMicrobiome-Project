#!/usr/bin/env python3
# Author: Lin
# Date: 2026-03-27
"""
Transfer learning comparison using Neural Network (MLPClassifier).

In each LODO fold, three strategies are evaluated on the same 25% test split
across 20 paired repeats:
  1. LODO:             train on N-1 → test on 25%
  2. Transfer:         pretrain on N-1 → fine-tune on 75% → test on 25%
  3. Transfer+ComBat:  same as Transfer but source goes through ComBat first

Per-fold signature genera are read from lodo_nested_*.py outputs (no leakage).
Architecture selected via inner 3-fold CV over shallow/medium/deep candidates.

outputs go to tables/:
  transfer_learning_comprehensive_results.csv
  transfer_architecture_results.csv
"""

import os
import copy
import numpy as np
import pandas as pd
from collections import Counter
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score, confusion_matrix
from pycombat import Combat
import warnings
warnings.filterwarnings("ignore")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = "/Users/lindan/Dropbox/PhD/Projects/Animal_Microbiome_cb"
TAB_DIR     = os.path.join(SCRIPT_DIR, "tables")
os.makedirs(TAB_DIR, exist_ok=True)

N_REPEATS    = 20
RANDOM_STATE = 42

# architectures and fixed hyperparameters
ARCHITECTURES = {
    "shallow": (32,),
    "medium":  (64, 32),
    "deep":    (64, 32, 16),
}
FIXED_ALPHA = 0.001
FIXED_LR    = 0.0001

# host configs
HOST_CONFIG = {
    "bovine": {
        "base_dir":           os.path.join(PROJECT_DIR, "metadata/bovine"),
        "studies":            ["PRJNA716761", "PRJNA1263132",
                               "PRJEB25741",  "PRJNA918869"],
        "min_study_presence": 3,
        "exclude_genera":     {"Dietzia", "Corynebacterium", "Roseburia",
                               "Anaerostipes", "Fournierella"},
    },
    "canine": {
        "base_dir":           os.path.join(PROJECT_DIR, "metadata/canine"),
        "studies":            ["PRJEB13362", "PRJNA905458", "PRJNA861113"],
        "min_study_presence": 2,
        "exclude_genera":     set(),
    },
}


def load_study_data(base_dir, study_id):
    """load otu + metadata, filter to normal/sick only"""
    study_dir = os.path.join(base_dir, study_id)
    otu  = pd.read_csv(os.path.join(study_dir, "otu_genus_update.csv"),
                       index_col=0)
    meta = pd.read_csv(os.path.join(study_dir, f"{study_id}.csv"))
    if "Run" not in meta.columns and "Unnamed: 0" in meta.columns:
        meta = meta.rename(columns={"Unnamed: 0": "Run"})
    meta   = meta[meta["group"].isin(["normal", "sick"])].copy()
    shared = list(set(otu.columns) & set(meta["Run"]))
    otu    = otu[shared]
    meta   = meta[meta["Run"].isin(shared)].copy()
    return otu, meta


def prepare_Xy_with_batch(study_data, study_ids, genera):
    """stack samples; applies per-study CZM (0.65×col-min) + CLR; returns batch labels for ComBat"""
    Xs, ys, batches = [], [], []
    for sid in study_ids:
        otu, meta = study_data[sid]
        samples = meta["Run"].values
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
                mat[:, j] = 1e-10
        # CLR per sample
        log_mat = np.log(mat)
        clr = log_mat - log_mat.mean(axis=1, keepdims=True)
        clr = np.nan_to_num(clr, nan=0.0, posinf=0.0, neginf=0.0)
        y = (meta.set_index("Run").loc[samples, "group"] == "sick"
             ).astype(int).values
        Xs.append(clr)
        ys.append(y)
        batches.extend([sid] * len(samples))
    return np.vstack(Xs), np.concatenate(ys), np.array(batches)


def get_fold_signatures(host):
    """load per-fold signature genera from lodo_nested_{host}.py output"""
    path = os.path.join(SCRIPT_DIR, host, "lodo_nested_fold_signatures.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Run lodo_nested_{host}.py first: {path}")
    df = pd.read_csv(path)
    return df.groupby("fold_test")["genus"].apply(sorted).to_dict()


def get_common_genera(study_data, studies, min_study_presence, exclude_genera):
    """common genera present in >= min_study_presence studies, minus exclusions"""
    genus_counts = Counter()
    for sid in studies:
        genus_counts.update(study_data[sid][0].index.tolist())
    return sorted([g for g, c in genus_counts.items()
                   if c >= min_study_presence and g not in exclude_genera])


# ComBat (source only, test untouched)
def apply_combat(X_train, batch_train):
    try:
        if len(np.unique(batch_train)) < 2:
            return X_train
        combat = Combat()
        X_c = combat.fit_transform(Y=X_train, b=batch_train)
        if np.any(~np.isfinite(X_c)):
            return X_train
        return X_c
    except Exception as e:
        print(f"    ComBat failed: {e}, using uncorrected data")
        return X_train


# NN helpers
def make_nn(arch, random_state=42):
    return MLPClassifier(
        hidden_layer_sizes=arch,
        alpha=FIXED_ALPHA,
        learning_rate_init=FIXED_LR,
        max_iter=500,
        early_stopping=True,
        validation_fraction=0.15,
        n_iter_no_change=20,
        random_state=random_state,
    )


def safe_metrics(y_true, y_prob):
    """AUC, sensitivity, specificity from probability scores; NaN-safe"""
    if len(np.unique(y_true)) < 2:
        return np.nan, np.nan, np.nan
    try:
        auroc = roc_auc_score(y_true, y_prob)
    except ValueError:
        auroc = np.nan
    y_pred = (y_prob >= 0.5).astype(int)
    try:
        tn, fp, fn, tp = confusion_matrix(
            y_true, y_pred, labels=[0, 1]).ravel()
        sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    except Exception:
        sens, spec = np.nan, np.nan
    return auroc, sens, spec


def cv_select_architecture(X_train, y_train, n_inner_folds=3, random_state=42):
    """pick best architecture via inner stratified CV on source training data"""
    n_min = min(np.bincount(y_train))
    n_splits = min(n_inner_folds, n_min)
    if n_splits < 2:
        return (64, 32), "medium"

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True,
                          random_state=random_state)
    best_score, best_arch, best_name = -1, (64, 32), "medium"

    for arch_name, arch in ARCHITECTURES.items():
        fold_aucs = []
        for tr_idx, val_idx in skf.split(X_train, y_train):
            nn   = make_nn(arch, random_state)
            nn.fit(X_train[tr_idx], y_train[tr_idx])
            prob = (nn.predict_proba(X_train[val_idx])[:, 1]
                    if hasattr(nn, "predict_proba")
                    else nn.decision_function(X_train[val_idx]))
            auc, _, _ = safe_metrics(y_train[val_idx], prob)
            if not np.isnan(auc):
                fold_aucs.append(auc)
        mean_auc = np.mean(fold_aucs) if fold_aucs else 0
        if mean_auc > best_score:
            best_score, best_arch, best_name = mean_auc, arch, arch_name

    return best_arch, best_name


# paired evaluation: all 3 strategies share the same 25% test set per repeat
def paired_eval(X_train, y_train, X_train_cb, X_target, y_target,
                arch, arch_cb):
    """evaluate all 3 strategies on 20 paired 75/25 splits"""
    sss = StratifiedShuffleSplit(n_splits=N_REPEATS, test_size=0.25,
                                  random_state=RANDOM_STATE)
    keys = ["lodo", "transfer", "transfer_combat"]
    out  = {k: {"auroc": [], "sensitivity": [], "specificity": []}
            for k in keys}

    for rep, (ft_idx, te_idx) in enumerate(sss.split(X_target, y_target)):
        if (len(np.unique(y_target[ft_idx])) < 2 or
                len(np.unique(y_target[te_idx])) < 2):
            continue

        X_ft, y_ft = X_target[ft_idx], y_target[ft_idx]
        X_te, y_te = X_target[te_idx], y_target[te_idx]
        rs = RANDOM_STATE + rep

        # LODO: train on source, test on same 25%
        nn = make_nn(arch, random_state=rs)
        nn.fit(X_train, y_train)
        prob = (nn.predict_proba(X_te)[:, 1]
                if hasattr(nn, "predict_proba")
                else nn.decision_function(X_te))
        auc, sens, spec = safe_metrics(y_te, prob)
        out["lodo"]["auroc"].append(auc)
        out["lodo"]["sensitivity"].append(sens)
        out["lodo"]["specificity"].append(spec)

        # Transfer (no ComBat): pretrain → fine-tune 75% → test 25%
        nn = make_nn(arch, random_state=rs)
        nn.fit(X_train, y_train)
        nn.set_params(warm_start=True, max_iter=200)
        nn.fit(X_ft, y_ft)
        prob = (nn.predict_proba(X_te)[:, 1]
                if hasattr(nn, "predict_proba")
                else nn.decision_function(X_te))
        auc, sens, spec = safe_metrics(y_te, prob)
        out["transfer"]["auroc"].append(auc)
        out["transfer"]["sensitivity"].append(sens)
        out["transfer"]["specificity"].append(spec)

        # Transfer + ComBat: pretrain on ComBat source → fine-tune → test
        nn = make_nn(arch_cb, random_state=rs)
        nn.fit(X_train_cb, y_train)
        nn.set_params(warm_start=True, max_iter=200)
        nn.fit(X_ft, y_ft)
        prob = (nn.predict_proba(X_te)[:, 1]
                if hasattr(nn, "predict_proba")
                else nn.decision_function(X_te))
        auc, sens, spec = safe_metrics(y_te, prob)
        out["transfer_combat"]["auroc"].append(auc)
        out["transfer_combat"]["sensitivity"].append(sens)
        out["transfer_combat"]["specificity"].append(spec)

    return out


# architecture comparison: transfer only, same 20 splits, no inner CV
def arch_eval(X_train, y_train, X_target, y_target):
    """evaluate all 3 architectures under transfer (no ComBat), signature features"""
    sss = StratifiedShuffleSplit(n_splits=N_REPEATS, test_size=0.25,
                                  random_state=RANDOM_STATE)
    rows = []
    for rep, (ft_idx, te_idx) in enumerate(sss.split(X_target, y_target)):
        if (len(np.unique(y_target[ft_idx])) < 2 or
                len(np.unique(y_target[te_idx])) < 2):
            continue
        X_ft, y_ft = X_target[ft_idx], y_target[ft_idx]
        X_te, y_te = X_target[te_idx], y_target[te_idx]
        rs = RANDOM_STATE + rep
        for arch_name, arch in ARCHITECTURES.items():
            nn = make_nn(arch, random_state=rs)
            nn.fit(X_train, y_train)
            nn.set_params(warm_start=True, max_iter=200)
            nn.fit(X_ft, y_ft)
            prob = (nn.predict_proba(X_te)[:, 1]
                    if hasattr(nn, "predict_proba")
                    else nn.decision_function(X_te))
            auc, sens, spec = safe_metrics(y_te, prob)
            rows.append(dict(architecture=arch_name, repeat=rep,
                             auroc=auc, sensitivity=sens, specificity=spec))
    return rows


def run_host(host):
    cfg      = HOST_CONFIG[host]
    studies  = cfg["studies"]
    base_dir = cfg["base_dir"]

    study_data = {}
    for sid in studies:
        otu, meta = load_study_data(base_dir, sid)
        study_data[sid] = (otu, meta)
        print(f"  {sid}: sick={sum(meta['group']=='sick')}, "
              f"normal={sum(meta['group']=='normal')}")

    # per-fold signatures
    fold_signatures = get_fold_signatures(host)
    print(f"  Fold signatures loaded for {len(fold_signatures)} folds")

    # common genera (full feature set)
    common_genera = get_common_genera(
        study_data, studies,
        cfg["min_study_presence"], cfg["exclude_genera"])
    print(f"  Common genera (full): {len(common_genera)}")

    results = []
    arch_results = []
    STRAT_META = {
        "lodo":            ("lodo",     False),
        "transfer":        ("transfer", False),
        "transfer_combat": ("transfer", True),
    }
    STRAT_LABEL = {
        "lodo": "LODO", "transfer": "Transfer",
        "transfer_combat": "Transfer+ComBat",
    }

    for test_study in studies:
        train_studies = [s for s in studies if s != test_study]

        for feature_set, genera in [
            ("signature", fold_signatures.get(test_study, [])),
            ("full",      common_genera),
        ]:
            if not genera:
                print(f"\n  --- {test_study} [{feature_set}]: "
                      f"no genera, skipping ---")
                continue

            print(f"\n  --- Test: {test_study}  [{feature_set}]  "
                  f"n_genera={len(genera)} ---")

            X_test,  y_test,  _ = prepare_Xy_with_batch(
                study_data, [test_study], genera)
            X_train, y_train, batch_train = prepare_Xy_with_batch(
                study_data, train_studies, genera)

            n_test_min = (min(np.bincount(y_test))
                          if len(np.unique(y_test)) >= 2 else 0)
            print(f"    Train: {len(y_train)} | "
                  f"Test: {len(y_test)} (min class: {n_test_min})")

            if n_test_min < 2:
                print(f"    SKIPPED (min class < 2)")
                continue

            X_train_cb = apply_combat(X_train, batch_train)

            best_arch,    arch_name    = cv_select_architecture(
                X_train,    y_train)
            best_arch_cb, arch_name_cb = cv_select_architecture(
                X_train_cb, y_train)
            print(f"    Arch (raw): {arch_name}  |  "
                  f"Arch (ComBat): {arch_name_cb}")

            paired = paired_eval(X_train, y_train, X_train_cb,
                                 X_test, y_test, best_arch, best_arch_cb)

            # arch comparison (signature only)
            if feature_set == "signature":
                for row in arch_eval(X_train, y_train, X_test, y_test):
                    row.update(host=host, test_study=test_study,
                               n_genera=len(genera))
                    arch_results.append(row)

            arch_map = {
                "lodo":            arch_name,
                "transfer":        arch_name,
                "transfer_combat": arch_name_cb,
            }
            for key in ["lodo", "transfer", "transfer_combat"]:
                strategy, combat = STRAT_META[key]
                metrics = paired[key]
                n = len(metrics["auroc"])
                mean_auc = np.nanmean(metrics["auroc"]) if n else np.nan
                print(f"    {STRAT_LABEL[key]:20s}: "
                      f"AUC={mean_auc:.3f} (n={n})")
                for i in range(n):
                    results.append(dict(
                        host=host,
                        test_study=test_study,
                        feature_set=feature_set,
                        strategy=strategy,
                        combat=combat,
                        best_architecture=arch_map[key],
                        n_genera=len(genera),
                        auroc=metrics["auroc"][i],
                        sensitivity=metrics["sensitivity"][i],
                        specificity=metrics["specificity"][i],
                    ))

    return pd.DataFrame(results), pd.DataFrame(arch_results)


if __name__ == "__main__":
    print("Transfer learning | bovine + canine")
    print(f"  strategies: LODO | Transfer | Transfer+ComBat")
    print(f"  architectures: shallow(32), medium(64,32), deep(64,32,16)")
    print(f"  fixed HP: alpha={FIXED_ALPHA}, lr={FIXED_LR}, repeats={N_REPEATS}")

    all_results  = []
    all_arch     = []
    for host in ["bovine", "canine"]:
        print(f"\n--- {host.upper()} ---")
        df, arch_df = run_host(host)
        all_results.append(df)
        all_arch.append(arch_df)

    results_df = pd.concat(all_results, ignore_index=True)
    out_path = os.path.join(TAB_DIR,
                            "transfer_learning_comprehensive_results.csv")
    results_df.to_csv(out_path, index=False)
    print(f"\nResults saved: {out_path}  ({len(results_df)} rows)")

    arch_df_all = pd.concat(all_arch, ignore_index=True)
    arch_path = os.path.join(TAB_DIR, "transfer_architecture_results.csv")
    arch_df_all.to_csv(arch_path, index=False)
    print(f"Architecture results saved: {arch_path}  ({len(arch_df_all)} rows)")

    print()
    for host in ["bovine", "canine"]:
        print(f"\n{host.upper()}:")
        df_h = results_df[results_df["host"] == host]
        for ts in df_h["test_study"].unique():
            df_ts = df_h[df_h["test_study"] == ts]
            print(f"  {ts}:")
            for fs in ["signature", "full"]:
                df_fs = df_ts[df_ts["feature_set"] == fs]
                for strat, cb in [("lodo", False), ("transfer", False),
                                   ("transfer", True)]:
                    sub = df_fs[(df_fs["strategy"] == strat) &
                                (df_fs["combat"] == cb)]["auroc"].dropna()
                    label = ("LODO" if strat == "lodo"
                             else ("Transfer+ComBat" if cb else "Transfer"))
                    if len(sub) > 0:
                        print(f"    [{fs:9s}] {label:20s}: "
                              f"{sub.mean()*100:.1f}+/-{sub.std()*100:.1f} "
                              f"(n={len(sub)})")
