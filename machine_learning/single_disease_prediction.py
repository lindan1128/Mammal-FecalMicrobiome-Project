#!/usr/bin/env python3
"""
Single Disease Prediction using Multiple Machine Learning Methods

This script performs disease prediction using three different methods:
- Random Forest
- Lasso Logistic Regression
- Neural Network

The script takes an OTU table and metadata as input, and performs
cross-validation for each method.
"""

import argparse
import sys
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score,
    accuracy_score,
    recall_score,
    f1_score,
    confusion_matrix,
    roc_curve
)
from sklearn.utils.class_weight import compute_class_weight
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader


def extract_genus(taxon_names):
    """
    Extract genus level classification from taxon names.
    
    Parameters:
    -----------
    taxon_names : list or pd.Index
        Taxon names to process
    
    Returns:
    --------
    list
        Extracted genus names
    """
    result = []
    
    for taxon_name in taxon_names:
        # Split by semicolon and get the 6th element (genus level)
        parts = taxon_name.split(';')
        
        if len(parts) >= 6:
            genus = parts[5]  # Index 5 is the 6th element (0-indexed)
            
            # Remove g__ prefix if present
            if genus.startswith('g__'):
                genus = genus[3:]  # Remove first 3 characters
            
            result.append(genus)
        else:
            # If not enough levels, return original
            result.append(taxon_name)
    
    return result


def rf_cv(data, labels, n_trees=20, n_folds=5, n_repeats=20, random_state=42):
    """
    Random forest cross-validation with class weights and feature importance.
    
    Parameters:
    -----------
    data : pd.DataFrame or np.ndarray
        Feature data (samples x features)
    labels : pd.Series or np.ndarray
        Class labels
    n_trees : int
        Number of trees in random forest
    n_folds : int
        Number of folds for cross-validation
    n_repeats : int
        Number of repeated cross-validation
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    dict with:
        - metrics_df: DataFrame with columns [metrics, value, repeats, folds]
        - avg_importance: Average feature importance across all repeats
        - feature_names: Names of features (if data is DataFrame)
    """
    print(f"Running random forest with {n_trees} trees, "
          f"{n_folds} folds, {n_repeats} repeats")
    print("Using balanced class weights to address class imbalance")
    
    # Convert to numpy arrays if needed
    if isinstance(data, pd.DataFrame):
        feature_names = data.columns.tolist()
        X = data.values
    else:
        feature_names = [f"feature_{i}" for i in range(data.shape[1])]
        X = data
    
    if isinstance(labels, pd.Series):
        y = labels.values
    else:
        y = labels
    
    # Check class distribution
    unique, counts = np.unique(y, return_counts=True)
    print("Original class distribution:")
    for label, count in zip(unique, counts):
        print(f"  {label}: {count}")
    
    # Initialize result storage
    n_total_folds = n_repeats * n_folds
    aucs_list = np.zeros(n_total_folds)
    acc_list = np.zeros(n_total_folds)
    recall_list = np.zeros(n_total_folds)
    f1_list = np.zeros(n_total_folds)
    
    importance_list = []
    
    fold_counter = 0
    
    for rep in range(n_repeats):
        print(f"\n=== Repeat {rep + 1} ===")
        
        # Create stratified k-fold
        skf = StratifiedKFold(
            n_splits=n_folds,
            shuffle=True,
            random_state=random_state + rep
        )
        
        fold_importance = []
        
        for fold, (train_idx, test_idx) in enumerate(skf.split(X, y)):
            print(f"  Fold {fold + 1}")
            
            # Split data
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            
            print(f"    Training samples: {len(train_idx)}, "
                  f"Test samples: {len(test_idx)}")
            train_unique, train_counts = np.unique(y_train, return_counts=True)
            test_unique, test_counts = np.unique(y_test, return_counts=True)
            print(f"    Training - {dict(zip(train_unique, train_counts))}")
            print(f"    Test - {dict(zip(test_unique, test_counts))}")
            
            # Train random forest with balanced class weights
            rf_model = RandomForestClassifier(
                n_estimators=n_trees,
                class_weight='balanced',
                random_state=random_state + rep,
                n_jobs=-1
            )
            
            rf_model.fit(X_train, y_train)
            
            # Predictions
            y_pred_proba = rf_model.predict_proba(X_test)[:, 1]
            y_pred = rf_model.predict(X_test)
            
            # Calculate metrics
            try:
                auc_val = roc_auc_score(y_test, y_pred_proba)
            except Exception:
                auc_val = np.nan
            
            acc_val = accuracy_score(y_test, y_pred)
            
            # For binary classification, use binary recall and f1
            recall_val = recall_score(
                y_test,
                y_pred,
                average='binary',
                pos_label=unique[1]
            )
            f1_val = f1_score(
                y_test,
                y_pred,
                average='binary',
                pos_label=unique[1]
            )
            
            # Store results
            aucs_list[fold_counter] = auc_val
            acc_list[fold_counter] = acc_val
            recall_list[fold_counter] = recall_val
            f1_list[fold_counter] = f1_val
            
            # Store feature importance
            fold_importance.append(rf_model.feature_importances_)
            
            fold_counter += 1
            
            print(f"    AUC: {auc_val:.3f}, Acc: {acc_val:.3f}, "
                  f"Recall: {recall_val:.3f}, F1: {f1_val:.3f}")
        
        # Average importance for this repeat
        repeat_avg_importance = np.mean(fold_importance, axis=0)
        importance_list.append(repeat_avg_importance)
        
        # Calculate mean AUC for this repeat
        repeat_start = rep * n_folds
        repeat_end = (rep + 1) * n_folds
        repeat_mean_auc = np.nanmean(aucs_list[repeat_start:repeat_end])
        print(f"  Repeat {rep + 1} completed - Mean AUC: {repeat_mean_auc:.3f}")
    
    # Calculate overall average importance
    avg_importance = np.mean(importance_list, axis=0)
    
    # Create metrics dataframe
    repeats_vec = np.repeat(np.arange(1, n_repeats + 1), n_folds)
    folds_vec = np.tile(np.arange(1, n_folds + 1), n_repeats)
    
    metrics_df = pd.DataFrame({
        'metrics': np.repeat(['auc', 'accuracy', 'recall', 'f1'],
                             n_total_folds),
        'value': np.concatenate([aucs_list, acc_list, recall_list, f1_list]),
        'repeats': np.tile(repeats_vec, 4),
        'folds': np.tile(folds_vec, 4)
    })
    
    print(f"\nCross-validation completed!")
    print(f"Metrics dataframe shape: {metrics_df.shape}")
    print("\nOverall summary:")
    summary = metrics_df.groupby('metrics')['value'].agg(['mean', 'std'])
    print(summary)
    
    # Create importance dataframe
    importance_df = pd.DataFrame({
        'feature': feature_names,
        'importance': avg_importance
    }).sort_values('importance', ascending=False)
    
    return {
        'metrics_df': metrics_df,
        'avg_importance': importance_df,
        'feature_names': feature_names
    }


def lasso_cv(data, labels, alpha=1.0, n_folds=5, n_repeats=20,
             random_state=42):
    """
    Lasso logistic regression cross-validation with class weights.
    
    Parameters:
    -----------
    data : pd.DataFrame or np.ndarray
        Feature data (samples x features)
    labels : pd.Series or np.ndarray
        Class labels
    alpha : float
        Regularization strength (C = 1/alpha in sklearn)
    n_folds : int
        Number of folds for cross-validation
    n_repeats : int
        Number of repeated cross-validation
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    dict with:
        - metrics_df: DataFrame with columns [metrics, value, repeats, folds]
    """
    print(f"Running Lasso logistic regression with alpha={alpha}, "
          f"{n_folds} folds, {n_repeats} repeats")
    print("Using balanced class weights to address class imbalance")
    
    # Convert to numpy arrays if needed
    if isinstance(data, pd.DataFrame):
        X = data.values
    else:
        X = data
    
    if isinstance(labels, pd.Series):
        y = labels.values
    else:
        y = labels
    
    # Check class distribution
    unique, counts = np.unique(y, return_counts=True)
    print("Original class distribution:")
    for label, count in zip(unique, counts):
        print(f"  {label}: {count}")
    
    # Initialize result storage
    n_total_folds = n_repeats * n_folds
    aucs_list = np.zeros(n_total_folds)
    acc_list = np.zeros(n_total_folds)
    recall_list = np.zeros(n_total_folds)
    f1_list = np.zeros(n_total_folds)
    
    fold_counter = 0
    
    for rep in range(n_repeats):
        print(f"\n=== Repeat {rep + 1} ===")
        
        # Create stratified k-fold
        skf = StratifiedKFold(
            n_splits=n_folds,
            shuffle=True,
            random_state=random_state + rep
        )
        
        for fold, (train_idx, test_idx) in enumerate(skf.split(X, y)):
            print(f"  Fold {fold + 1}")
            
            # Split data
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            
            print(f"    Training samples: {len(train_idx)}, "
                  f"Test samples: {len(test_idx)}")
            train_unique, train_counts = np.unique(y_train, return_counts=True)
            test_unique, test_counts = np.unique(y_test, return_counts=True)
            print(f"    Training - {dict(zip(train_unique, train_counts))}")
            print(f"    Test - {dict(zip(test_unique, test_counts))}")
            
            # Standardize features (important for Lasso)
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)
            
            # Train Lasso logistic regression with balanced class weights
            lasso_model = LogisticRegression(
                penalty='l1',
                C=1.0 / alpha,
                class_weight='balanced',
                solver='liblinear',
                random_state=random_state + rep,
                max_iter=1000
            )
            
            lasso_model.fit(X_train_scaled, y_train)
            
            # Predictions
            y_pred_proba = lasso_model.predict_proba(X_test_scaled)[:, 1]
            y_pred = lasso_model.predict(X_test_scaled)
            
            # Calculate metrics
            try:
                auc_val = roc_auc_score(y_test, y_pred_proba)
            except Exception:
                auc_val = np.nan
            
            acc_val = accuracy_score(y_test, y_pred)
            
            # For binary classification, use binary recall and f1
            recall_val = recall_score(
                y_test,
                y_pred,
                average='binary',
                pos_label=unique[1]
            )
            f1_val = f1_score(
                y_test,
                y_pred,
                average='binary',
                pos_label=unique[1]
            )
            
            # Store results
            aucs_list[fold_counter] = auc_val
            acc_list[fold_counter] = acc_val
            recall_list[fold_counter] = recall_val
            f1_list[fold_counter] = f1_val
            
            fold_counter += 1
            
            print(f"    AUC: {auc_val:.3f}, Acc: {acc_val:.3f}, "
                  f"Recall: {recall_val:.3f}, F1: {f1_val:.3f}")
        
        # Calculate mean AUC for this repeat
        repeat_start = rep * n_folds
        repeat_end = (rep + 1) * n_folds
        repeat_mean_auc = np.nanmean(aucs_list[repeat_start:repeat_end])
        print(f"  Repeat {rep + 1} completed - Mean AUC: {repeat_mean_auc:.3f}")
    
    # Create metrics dataframe
    repeats_vec = np.repeat(np.arange(1, n_repeats + 1), n_folds)
    folds_vec = np.tile(np.arange(1, n_folds + 1), n_repeats)
    
    metrics_df = pd.DataFrame({
        'metrics': np.repeat(['auc', 'accuracy', 'recall', 'f1'],
                             n_total_folds),
        'value': np.concatenate([aucs_list, acc_list, recall_list, f1_list]),
        'repeats': np.tile(repeats_vec, 4),
        'folds': np.tile(folds_vec, 4)
    })
    
    print(f"\nCross-validation completed!")
    print(f"Metrics dataframe shape: {metrics_df.shape}")
    print("\nOverall summary:")
    summary = metrics_df.groupby('metrics')['value'].agg(['mean', 'std'])
    print(summary)
    
    return {
        'metrics_df': metrics_df
    }


class SimpleNN(nn.Module):
    """Simple Neural Network for binary classification."""
    
    def __init__(self, input_size):
        super(SimpleNN, self).__init__()
        self.layer1 = nn.Linear(input_size, 16)
        self.layer2 = nn.Linear(16, 8)
        self.layer3 = nn.Linear(8, 1)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        
    def forward(self, x):
        x = self.relu(self.layer1(x))
        x = self.relu(self.layer2(x))
        x = self.sigmoid(self.layer3(x))
        return x


def nn_cv(data, labels, n_epochs=100, batch_size=32, lr=0.001,
          n_folds=5, n_repeats=20, random_state=42):
    """
    Neural network cross-validation with class weights.
    
    Parameters:
    -----------
    data : pd.DataFrame or np.ndarray
        Feature data (samples x features)
    labels : pd.Series or np.ndarray
        Class labels
    n_epochs : int
        Number of training epochs
    batch_size : int
        Batch size for training
    lr : float
        Learning rate
    n_folds : int
        Number of folds for cross-validation
    n_repeats : int
        Number of repeated cross-validation
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    dict with:
        - metrics_df: DataFrame with columns [metrics, value, repeats, folds]
    """
    print(f"Running Neural Network with {n_epochs} epochs, "
          f"{n_folds} folds, {n_repeats} repeats")
    print(f"Batch size: {batch_size}, Learning rate: {lr}")
    print("Using weighted loss to address class imbalance")
    
    # Set random seeds for reproducibility
    torch.manual_seed(random_state)
    np.random.seed(random_state)
    
    # Convert to numpy arrays if needed
    if isinstance(data, pd.DataFrame):
        X = data.values
    else:
        X = data
    
    if isinstance(labels, pd.Series):
        y = labels.values
    else:
        y = labels
    
    # Check class distribution
    unique, counts = np.unique(y, return_counts=True)
    print("Original class distribution:")
    for label, count in zip(unique, counts):
        print(f"  {label}: {count}")
    
    # Convert labels to binary (0, 1)
    label_mapping = {unique[0]: 0, unique[1]: 1}
    y_binary = np.array([label_mapping[label] for label in y])
    
    # Calculate class weights
    total_samples = len(y_binary)
    class_counts = np.bincount(y_binary)
    class_weights = total_samples / (len(class_counts) * class_counts)
    print(f"Class weights: {dict(zip([0, 1], class_weights))}")
    
    # Initialize result storage
    n_total_folds = n_repeats * n_folds
    aucs_list = np.zeros(n_total_folds)
    acc_list = np.zeros(n_total_folds)
    recall_list = np.zeros(n_total_folds)
    f1_list = np.zeros(n_total_folds)
    
    fold_counter = 0
    
    # Set device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    for rep in range(n_repeats):
        print(f"\n=== Repeat {rep + 1} ===")
        
        # Create stratified k-fold
        skf = StratifiedKFold(
            n_splits=n_folds,
            shuffle=True,
            random_state=random_state + rep
        )
        
        for fold, (train_idx, test_idx) in enumerate(skf.split(X, y_binary)):
            print(f"  Fold {fold + 1}")
            
            # Split data
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y_binary[train_idx], y_binary[test_idx]
            
            print(f"    Training samples: {len(train_idx)}, "
                  f"Test samples: {len(test_idx)}")
            train_unique, train_counts = np.unique(y_train, return_counts=True)
            test_unique, test_counts = np.unique(y_test, return_counts=True)
            print(f"    Training - {dict(zip(train_unique, train_counts))}")
            print(f"    Test - {dict(zip(test_unique, test_counts))}")
            
            # Standardize features
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)
            
            # Convert to tensors
            X_train_tensor = torch.FloatTensor(X_train_scaled).to(device)
            y_train_tensor = torch.FloatTensor(y_train).unsqueeze(1).to(device)
            X_test_tensor = torch.FloatTensor(X_test_scaled).to(device)
            
            # Create dataset and dataloader
            train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
            train_loader = DataLoader(
                train_dataset,
                batch_size=batch_size,
                shuffle=True
            )
            
            # Initialize model
            input_size = X_train.shape[1]
            model = SimpleNN(input_size).to(device)
            
            # Calculate sample weights for weighted loss
            sample_weights = torch.FloatTensor(
                [class_weights[int(label)] for label in y_train]
            ).to(device)
            
            # Loss function and optimizer
            criterion = nn.BCELoss(reduction='none')
            optimizer = optim.Adam(model.parameters(), lr=lr)
            
            # Training
            model.train()
            for epoch in range(n_epochs):
                epoch_loss = 0.0
                for batch_X, batch_y in train_loader:
                    # Forward pass
                    outputs = model(batch_X)
                    
                    # Calculate weighted loss
                    batch_indices = []
                    for i in range(len(batch_y)):
                        idx = (
                            (X_train_tensor == batch_X[i]).all(dim=1)
                        ).nonzero(as_tuple=True)[0]
                        if len(idx) > 0:
                            batch_indices.append(idx[0].item())
                    
                    if len(batch_indices) > 0:
                        batch_weights = sample_weights[batch_indices]
                        loss = criterion(outputs, batch_y)
                        weighted_loss = (
                            loss * batch_weights.unsqueeze(1)
                        ).mean()
                    else:
                        loss = criterion(outputs, batch_y).mean()
                        weighted_loss = loss
                    
                    # Backward pass
                    optimizer.zero_grad()
                    weighted_loss.backward()
                    optimizer.step()
                    
                    epoch_loss += weighted_loss.item()
            
            # Evaluation
            model.eval()
            with torch.no_grad():
                y_pred_proba = model(X_test_tensor).cpu().numpy().flatten()
                y_pred = (y_pred_proba >= 0.5).astype(int)
            
            # Calculate metrics
            try:
                auc_val = roc_auc_score(y_test, y_pred_proba)
            except Exception:
                auc_val = np.nan
            
            acc_val = accuracy_score(y_test, y_pred)
            recall_val = recall_score(
                y_test,
                y_pred,
                average='binary',
                pos_label=1
            )
            f1_val = f1_score(y_test, y_pred, average='binary', pos_label=1)
            
            # Store results
            aucs_list[fold_counter] = auc_val
            acc_list[fold_counter] = acc_val
            recall_list[fold_counter] = recall_val
            f1_list[fold_counter] = f1_val
            
            fold_counter += 1
            
            print(f"    AUC: {auc_val:.3f}, Acc: {acc_val:.3f}, "
                  f"Recall: {recall_val:.3f}, F1: {f1_val:.3f}")
        
        # Calculate mean AUC for this repeat
        repeat_start = rep * n_folds
        repeat_end = (rep + 1) * n_folds
        repeat_mean_auc = np.nanmean(aucs_list[repeat_start:repeat_end])
        print(f"  Repeat {rep + 1} completed - Mean AUC: {repeat_mean_auc:.3f}")
    
    # Create metrics dataframe
    repeats_vec = np.repeat(np.arange(1, n_repeats + 1), n_folds)
    folds_vec = np.tile(np.arange(1, n_folds + 1), n_repeats)
    
    metrics_df = pd.DataFrame({
        'metrics': np.repeat(['auc', 'accuracy', 'recall', 'f1'],
                             n_total_folds),
        'value': np.concatenate([aucs_list, acc_list, recall_list, f1_list]),
        'repeats': np.tile(repeats_vec, 4),
        'folds': np.tile(folds_vec, 4)
    })
    
    print(f"\nCross-validation completed!")
    print(f"Metrics dataframe shape: {metrics_df.shape}")
    print("\nOverall summary:")
    summary = metrics_df.groupby('metrics')['value'].agg(['mean', 'std'])
    print(summary)
    
    return {
        'metrics_df': metrics_df
    }


def load_and_preprocess_data(otu_file, metadata_file, group_column='group',
                              sample_column='Sample', target_groups=None,
                              abundance_threshold=0.001,
                              signature_set=None):
    """
    Load and preprocess OTU table and metadata.
    
    Parameters:
    -----------
    otu_file : str
        Path to OTU table CSV file
    metadata_file : str
        Path to metadata CSV file
    group_column : str
        Column name in metadata containing group labels
    sample_column : str
        Column name in metadata containing sample IDs
    target_groups : list or None
        List of target groups to include (if None, include all)
    abundance_threshold : float
        Minimum abundance threshold for filtering taxa
    signature_set : list or None
        List of signature taxa to use (if None, use all filtered taxa)
    
    Returns:
    --------
    dict with:
        - full_data: Full filtered OTU data
        - signature_data: Signature OTU data (if signature_set provided)
        - labels: Group labels
        - metadata: Processed metadata
    """
    print(f"Loading OTU data from: {otu_file}")
    otu_data = pd.read_csv(otu_file, index_col=0)
    
    print(f"Loading metadata from: {metadata_file}")
    metadata = pd.read_csv(metadata_file, index_col=0)
    
    # Filter metadata by target groups if specified
    if target_groups is not None:
        print(f"Filtering metadata for groups: {target_groups}")
        metadata = metadata[
            metadata[group_column].isin(target_groups)
        ].copy()
    
    # Find common samples
    common_samples = list(set(otu_data.columns) & set(metadata[sample_column]))
    print(f"Found {len(common_samples)} common samples")
    
    if len(common_samples) == 0:
        raise ValueError("No common samples found between OTU data and metadata")
    
    # Filter metadata and OTU data to common samples
    metadata = metadata[
        metadata[sample_column].isin(common_samples)
    ].copy()
    otu_data = otu_data[common_samples]
    
    # Transpose OTU data so rows=samples, columns=taxa
    otu_transposed = otu_data.T
    
    # Extract genus level
    print("Extracting genus level...")
    genus_names = extract_genus(otu_transposed.columns.tolist())
    otu_genus = otu_transposed.copy()
    otu_genus.columns = genus_names
    
    print(f"OTU data shape: {otu_genus.shape}")
    
    # Filter low abundance taxa
    print(f"Filtering taxa with abundance < {abundance_threshold}")
    col_sums = otu_genus.sum(axis=0)
    otu_filtered = otu_genus.loc[:, col_sums >= abundance_threshold]
    
    # Remove empty and uncultured taxa
    columns_to_keep = [
        col for col in otu_filtered.columns
        if col != '__' and 'uncultured' not in col.lower()
    ]
    otu_filtered = otu_filtered[columns_to_keep]
    
    print(f"After filtering: {otu_filtered.shape[1]} taxa remaining")
    
    # Prepare labels
    metadata = metadata.set_index(sample_column)
    metadata = metadata.loc[otu_filtered.index]
    
    # Get unique groups for categorical ordering
    unique_groups = metadata[group_column].unique().tolist()
    
    # Create categorical labels
    metadata['group_factor'] = pd.Categorical(
        metadata[group_column],
        categories=unique_groups
    )
    
    print(f"\nClass distribution:")
    print(metadata['group_factor'].value_counts())
    
    result = {
        'full_data': otu_filtered,
        'labels': metadata['group_factor'],
        'metadata': metadata
    }
    
    # Process signature data if provided
    if signature_set is not None:
        available_signature = [
            taxa for taxa in signature_set
            if taxa in otu_filtered.columns
        ]
        print(f"\nSignature taxa available in data: "
              f"{len(available_signature)}/{len(signature_set)}")
        
        if len(available_signature) > 0:
            signature_data = otu_filtered[available_signature]
            result['signature_data'] = signature_data
            print(f"Signature data shape: {signature_data.shape}")
    
    return result


def main():
    """Main function to run disease prediction analysis."""
    parser = argparse.ArgumentParser(
        description='Single Disease Prediction using Multiple ML Methods'
    )
    
    # Input files
    parser.add_argument(
        '--otu_file',
        type=str,
        required=True,
        help='Path to OTU table CSV file'
    )
    parser.add_argument(
        '--metadata_file',
        type=str,
        required=True,
        help='Path to metadata CSV file'
    )
    
    # Column specifications
    parser.add_argument(
        '--group_column',
        type=str,
        default='group',
        help='Column name in metadata containing group labels (default: group)'
    )
    parser.add_argument(
        '--sample_column',
        type=str,
        default='Sample',
        help='Column name in metadata containing sample IDs (default: Sample)'
    )
    
    # Target groups
    parser.add_argument(
        '--target_groups',
        type=str,
        nargs='+',
        default=None,
        help='Target groups to include in analysis (default: all groups)'
    )
    
    # Filtering parameters
    parser.add_argument(
        '--abundance_threshold',
        type=float,
        default=0.001,
        help='Minimum abundance threshold for filtering taxa (default: 0.001)'
    )
    
    # Signature set
    parser.add_argument(
        '--signature_set',
        type=str,
        nargs='+',
        default=None,
        help='List of signature taxa to use (default: None)'
    )
    
    # Model selection
    parser.add_argument(
        '--methods',
        type=str,
        nargs='+',
        choices=['rf', 'lasso', 'nn', 'all'],
        default=['all'],
        help='Methods to run: rf, lasso, nn, or all (default: all)'
    )
    
    # Cross-validation parameters
    parser.add_argument(
        '--n_folds',
        type=int,
        default=5,
        help='Number of folds for cross-validation (default: 5)'
    )
    parser.add_argument(
        '--n_repeats',
        type=int,
        default=20,
        help='Number of repeated cross-validation (default: 20)'
    )
    parser.add_argument(
        '--random_state',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    
    # Random Forest parameters
    parser.add_argument(
        '--n_trees',
        type=int,
        default=20,
        help='Number of trees in random forest (default: 20)'
    )
    
    # Lasso parameters
    parser.add_argument(
        '--alpha',
        type=float,
        default=1.0,
        help='Regularization strength for Lasso (default: 1.0)'
    )
    
    # Neural Network parameters
    parser.add_argument(
        '--n_epochs',
        type=int,
        default=100,
        help='Number of training epochs for neural network (default: 100)'
    )
    parser.add_argument(
        '--batch_size',
        type=int,
        default=32,
        help='Batch size for neural network training (default: 32)'
    )
    parser.add_argument(
        '--learning_rate',
        type=float,
        default=0.001,
        help='Learning rate for neural network (default: 0.001)'
    )
    
    # Output options
    parser.add_argument(
        '--output_prefix',
        type=str,
        default='results',
        help='Prefix for output files (default: results)'
    )
    parser.add_argument(
        '--save_results',
        action='store_true',
        help='Save results to CSV files'
    )
    
    args = parser.parse_args()
    
    # Determine which methods to run
    if 'all' in args.methods:
        methods_to_run = ['rf', 'lasso', 'nn']
    else:
        methods_to_run = args.methods
    
    print("="*70)
    print("Single Disease Prediction Analysis")
    print("="*70)
    print(f"Methods to run: {', '.join(methods_to_run)}")
    print(f"Cross-validation: {args.n_repeats} repeats x {args.n_folds} folds")
    print("="*70)
    
    # Load and preprocess data
    print("\n" + "="*70)
    print("LOADING AND PREPROCESSING DATA")
    print("="*70)
    
    data_dict = load_and_preprocess_data(
        otu_file=args.otu_file,
        metadata_file=args.metadata_file,
        group_column=args.group_column,
        sample_column=args.sample_column,
        target_groups=args.target_groups,
        abundance_threshold=args.abundance_threshold,
        signature_set=args.signature_set
    )
    
    full_data = data_dict['full_data']
    labels = data_dict['labels']
    signature_data = data_dict.get('signature_data', None)
    
    results = {}
    
    # Run Random Forest
    if 'rf' in methods_to_run:
        print("\n" + "="*70)
        print("RANDOM FOREST - FULL DATASET")
        print("="*70)
        
        rf_full_results = rf_cv(
            full_data,
            labels,
            n_trees=args.n_trees,
            n_folds=args.n_folds,
            n_repeats=args.n_repeats,
            random_state=args.random_state
        )
        results['rf_full'] = rf_full_results
        
        if args.save_results:
            rf_full_results['metrics_df'].to_csv(
                f"{args.output_prefix}_rf_full_metrics.csv",
                index=False
            )
            rf_full_results['avg_importance'].to_csv(
                f"{args.output_prefix}_rf_full_importance.csv",
                index=False
            )
            print(f"\nSaved Random Forest results to "
                  f"{args.output_prefix}_rf_full_*.csv")
        
        # Run on signature data if available
        if signature_data is not None:
            print("\n" + "="*70)
            print("RANDOM FOREST - SIGNATURE DATASET")
            print("="*70)
            
            rf_sig_results = rf_cv(
                signature_data,
                labels,
                n_trees=args.n_trees,
                n_folds=args.n_folds,
                n_repeats=args.n_repeats,
                random_state=args.random_state
            )
            results['rf_signature'] = rf_sig_results
            
            if args.save_results:
                rf_sig_results['metrics_df'].to_csv(
                    f"{args.output_prefix}_rf_signature_metrics.csv",
                    index=False
                )
                rf_sig_results['avg_importance'].to_csv(
                    f"{args.output_prefix}_rf_signature_importance.csv",
                    index=False
                )
                print(f"\nSaved Random Forest signature results to "
                      f"{args.output_prefix}_rf_signature_*.csv")
    
    # Run Lasso
    if 'lasso' in methods_to_run:
        print("\n" + "="*70)
        print("LASSO REGRESSION - FULL DATASET")
        print("="*70)
        
        lasso_full_results = lasso_cv(
            full_data,
            labels,
            alpha=args.alpha,
            n_folds=args.n_folds,
            n_repeats=args.n_repeats,
            random_state=args.random_state
        )
        results['lasso_full'] = lasso_full_results
        
        if args.save_results:
            lasso_full_results['metrics_df'].to_csv(
                f"{args.output_prefix}_lasso_full_metrics.csv",
                index=False
            )
            print(f"\nSaved Lasso results to "
                  f"{args.output_prefix}_lasso_full_metrics.csv")
        
        # Run on signature data if available
        if signature_data is not None:
            print("\n" + "="*70)
            print("LASSO REGRESSION - SIGNATURE DATASET")
            print("="*70)
            
            lasso_sig_results = lasso_cv(
                signature_data,
                labels,
                alpha=args.alpha,
                n_folds=args.n_folds,
                n_repeats=args.n_repeats,
                random_state=args.random_state
            )
            results['lasso_signature'] = lasso_sig_results
            
            if args.save_results:
                lasso_sig_results['metrics_df'].to_csv(
                    f"{args.output_prefix}_lasso_signature_metrics.csv",
                    index=False
                )
                print(f"\nSaved Lasso signature results to "
                      f"{args.output_prefix}_lasso_signature_metrics.csv")
    
    # Run Neural Network
    if 'nn' in methods_to_run:
        print("\n" + "="*70)
        print("NEURAL NETWORK - FULL DATASET")
        print("="*70)
        
        nn_full_results = nn_cv(
            full_data,
            labels,
            n_epochs=args.n_epochs,
            batch_size=args.batch_size,
            lr=args.learning_rate,
            n_folds=args.n_folds,
            n_repeats=args.n_repeats,
            random_state=args.random_state
        )
        results['nn_full'] = nn_full_results
        
        if args.save_results:
            nn_full_results['metrics_df'].to_csv(
                f"{args.output_prefix}_nn_full_metrics.csv",
                index=False
            )
            print(f"\nSaved Neural Network results to "
                  f"{args.output_prefix}_nn_full_metrics.csv")
        
        # Run on signature data if available
        if signature_data is not None:
            print("\n" + "="*70)
            print("NEURAL NETWORK - SIGNATURE DATASET")
            print("="*70)
            
            nn_sig_results = nn_cv(
                signature_data,
                labels,
                n_epochs=args.n_epochs,
                batch_size=args.batch_size,
                lr=args.learning_rate,
                n_folds=args.n_folds,
                n_repeats=args.n_repeats,
                random_state=args.random_state
            )
            results['nn_signature'] = nn_sig_results
            
            if args.save_results:
                nn_sig_results['metrics_df'].to_csv(
                    f"{args.output_prefix}_nn_signature_metrics.csv",
                    index=False
                )
                print(f"\nSaved Neural Network signature results to "
                      f"{args.output_prefix}_nn_signature_metrics.csv")
    
    # Print final summary
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70)
    
    for key, result in results.items():
        print(f"\n{key.upper().replace('_', ' ')}:")
        summary = result['metrics_df'].groupby('metrics')['value'].agg(
            ['mean', 'std']
        )
        print(summary)
    
    print("\n" + "="*70)
    print("Analysis completed successfully!")
    print("="*70)
    
    return results


if __name__ == "__main__":
    main()

