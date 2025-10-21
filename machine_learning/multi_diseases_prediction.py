#!/usr/bin/env python3
"""
Multi-Disease Prediction using Multiple Machine Learning Methods

This script performs multi-class disease prediction using:
- Random Forest
- Neural Network

The script supports batch effect correction using pyComBat and uses
One-vs-Rest strategy for per-class metrics.
"""

import argparse
import sys
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (
    roc_auc_score,
    accuracy_score,
    recall_score,
    f1_score,
    average_precision_score,
    confusion_matrix
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


def rf_cv_multiclass(data, labels, batch=None, apply_combat=False,
                      n_trees=20, n_folds=5, n_repeats=20, random_state=42):
    """
    Random forest cross-validation for multi-class classification.
    Uses One-vs-Rest strategy for per-class metrics including AUPR.
    
    Parameters:
    -----------
    data : pd.DataFrame or np.ndarray
        Feature data (samples x features)
    labels : pd.Series or np.ndarray
        Class labels (can be 'normal', 'diarrhea', or other diseases)
    batch : pd.Series or np.ndarray or None
        Batch labels for each sample (e.g., 'study1', 'study2')
        If None or apply_combat=False, no batch correction is applied
    apply_combat : bool
        Whether to apply batch effect correction using pyComBat
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
        - metrics_df: DataFrame with columns [disease, metrics, value,
                      repeats, folds]
        - avg_importance: Average feature importance across all repeats
        - feature_names: Names of features
    """
    print(f"Running random forest (multi-class) with {n_trees} trees, "
          f"{n_folds} folds, {n_repeats} repeats")
    print("Using balanced class weights to address class imbalance")
    
    if apply_combat and batch is not None:
        print("Batch effect correction: Using pyComBat")
        from pycombat import Combat
    else:
        print("Batch effect correction: None")
    
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
        y = np.array(labels)
    
    if batch is not None:
        if isinstance(batch, pd.Series):
            batch_array = batch.values
        else:
            batch_array = np.array(batch)
    else:
        batch_array = None
    
    # Get all classes
    classes = np.unique(y)
    n_classes = len(classes)
    print(f"Number of classes: {n_classes}")
    print("Class distribution:")
    for label in classes:
        count = np.sum(y == label)
        print(f"  {label}: {count}")
    
    if batch_array is not None:
        print("Batch distribution:")
        unique_batches, batch_counts = np.unique(
            batch_array,
            return_counts=True
        )
        for batch_name, count in zip(unique_batches, batch_counts):
            print(f"  {batch_name}: {count}")
    
    # Initialize result storage
    n_total_folds = n_repeats * n_folds
    
    # Store metrics for each class
    results = {cls: {
        'auc': np.zeros(n_total_folds),
        'aupr': np.zeros(n_total_folds),
        'accuracy': np.zeros(n_total_folds),
        'recall': np.zeros(n_total_folds),
        'f1': np.zeros(n_total_folds)
    } for cls in classes}
    
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
            
            # Apply batch effect correction if requested
            if apply_combat and batch_array is not None:
                train_batch = batch_array[train_idx]
                test_batch = batch_array[test_idx]
                
                # Check if we have multiple batches in training set
                unique_train_batches = np.unique(train_batch)
                
                if len(unique_train_batches) >= 2:
                    # Check if test batches are subset of train batches
                    unique_test_batches = np.unique(test_batch)
                    test_in_train = all(
                        b in unique_train_batches for b in unique_test_batches
                    )
                    
                    if test_in_train:
                        print(f"    Applying Combat correction "
                              f"({len(unique_train_batches)} batches)...")
                        
                        try:
                            # pyComBat expects: Y is (observations, features),
                            # b is (observations,)
                            # X_train is already (samples, features)
                            
                            # Initialize Combat
                            combat = Combat()
                            
                            # Fit on training data
                            X_train_corrected = combat.fit_transform(
                                Y=X_train,
                                b=train_batch
                            )
                            
                            # Transform test data
                            X_test_corrected = combat.transform(
                                Y=X_test,
                                b=test_batch
                            )
                            
                            # Check for NaN or Inf
                            if (np.isnan(X_train_corrected).any() or
                                    np.isinf(X_train_corrected).any()):
                                print("    Warning: Combat produced NaN/Inf, "
                                      "using original data")
                            else:
                                X_train = X_train_corrected
                                X_test = X_test_corrected
                                print("    Combat correction successful")
                        
                        except Exception as e:
                            print(f"    Error in ComBat: {str(e)}, "
                                  f"using original data")
                    else:
                        print("    Skipping ComBat: test batches not in "
                              "training set")
                else:
                    print(f"    Skipping ComBat: only "
                          f"{len(unique_train_batches)} batch in training set")
            
            print(f"    Training samples: {len(train_idx)}, "
                  f"Test samples: {len(test_idx)}")
            
            # Train random forest with balanced class weights
            rf_model = RandomForestClassifier(
                n_estimators=n_trees,
                class_weight='balanced',
                random_state=random_state + rep,
                n_jobs=-1
            )
            
            rf_model.fit(X_train, y_train)
            
            # Predictions
            y_pred = rf_model.predict(X_test)
            y_pred_proba = rf_model.predict_proba(X_test)
            
            # Overall accuracy (same for all classes)
            acc_val = accuracy_score(y_test, y_pred)
            
            # Calculate per-class metrics using One-vs-Rest
            for cls in classes:
                # Binary labels for this class (OvR)
                y_test_binary = (y_test == cls).astype(int)
                y_pred_binary = (y_pred == cls).astype(int)
                
                # Get probability for this class
                cls_idx = np.where(rf_model.classes_ == cls)[0][0]
                y_pred_proba_cls = y_pred_proba[:, cls_idx]
                
                # AUC (One-vs-Rest)
                try:
                    auc_val = roc_auc_score(y_test_binary, y_pred_proba_cls)
                except Exception:
                    auc_val = np.nan
                
                # AUPR (One-vs-Rest)
                try:
                    aupr_val = average_precision_score(
                        y_test_binary,
                        y_pred_proba_cls
                    )
                except Exception:
                    aupr_val = np.nan
                
                # Recall and F1 for this class
                recall_val = recall_score(
                    y_test_binary,
                    y_pred_binary,
                    zero_division=0
                )
                f1_val = f1_score(
                    y_test_binary,
                    y_pred_binary,
                    zero_division=0
                )
                
                # Store results
                results[cls]['auc'][fold_counter] = auc_val
                results[cls]['aupr'][fold_counter] = aupr_val
                results[cls]['accuracy'][fold_counter] = acc_val
                results[cls]['recall'][fold_counter] = recall_val
                results[cls]['f1'][fold_counter] = f1_val
            
            # Store feature importance
            fold_importance.append(rf_model.feature_importances_)
            fold_counter += 1
            
            print(f"    Overall Acc: {acc_val:.3f}")
            for cls in classes:
                print(f"    {cls} - "
                      f"AUC: {results[cls]['auc'][fold_counter-1]:.3f}, "
                      f"AUPR: {results[cls]['aupr'][fold_counter-1]:.3f}, "
                      f"Recall: {results[cls]['recall'][fold_counter-1]:.3f}, "
                      f"F1: {results[cls]['f1'][fold_counter-1]:.3f}")
        
        # Average importance for this repeat
        repeat_avg_importance = np.mean(fold_importance, axis=0)
        importance_list.append(repeat_avg_importance)
    
    # Calculate overall average importance
    avg_importance = np.mean(importance_list, axis=0)
    
    # Create metrics dataframe
    repeats_vec = np.repeat(np.arange(1, n_repeats + 1), n_folds)
    folds_vec = np.tile(np.arange(1, n_folds + 1), n_repeats)
    
    metrics_data = []
    for cls in classes:
        for metric in ['auc', 'aupr', 'accuracy', 'recall', 'f1']:
            for i in range(n_total_folds):
                metrics_data.append({
                    'disease': cls,
                    'metrics': metric,
                    'value': results[cls][metric][i],
                    'repeats': repeats_vec[i],
                    'folds': folds_vec[i]
                })
    
    metrics_df = pd.DataFrame(metrics_data)
    
    print(f"\nCross-validation completed!")
    print(f"Metrics dataframe shape: {metrics_df.shape}")
    print("\nOverall summary (per class):")
    summary = metrics_df.groupby(['disease', 'metrics'])['value'].agg(
        ['mean', 'std']
    )
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


class MultiClassNN(nn.Module):
    """Neural Network for multi-class classification."""
    
    def __init__(self, input_size, n_classes):
        super(MultiClassNN, self).__init__()
        self.layer1 = nn.Linear(input_size, 16)
        self.layer2 = nn.Linear(16, 8)
        self.layer3 = nn.Linear(8, n_classes)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(0.2)
        
    def forward(self, x):
        x = self.relu(self.layer1(x))
        x = self.dropout(x)
        x = self.relu(self.layer2(x))
        x = self.dropout(x)
        x = self.layer3(x)  # No activation (CrossEntropyLoss includes softmax)
        return x


def nn_cv_multiclass(data, labels, batch=None, apply_combat=False,
                      n_epochs=100, batch_size=32, lr=0.001,
                      n_folds=5, n_repeats=20, random_state=42):
    """
    Neural network cross-validation for multi-class classification.
    Uses One-vs-Rest strategy for per-class metrics including AUPR.
    
    Parameters:
    -----------
    data : pd.DataFrame or np.ndarray
        Feature data (samples x features)
    labels : pd.Series or np.ndarray
        Class labels (can be 'normal', 'diarrhea', or other diseases)
    batch : pd.Series or np.ndarray or None
        Batch labels for each sample (e.g., 'study1', 'study2')
        If None or apply_combat=False, no batch correction is applied
    apply_combat : bool
        Whether to apply batch effect correction using pyComBat
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
        - metrics_df: DataFrame with columns [disease, metrics, value,
                      repeats, folds]
    """
    print(f"Running Neural Network (multi-class) with {n_epochs} epochs, "
          f"{n_folds} folds, {n_repeats} repeats")
    print(f"Batch size: {batch_size}, Learning rate: {lr}")
    print("Using weighted loss to address class imbalance")
    
    if apply_combat and batch is not None:
        print("Batch effect correction: Using pyComBat")
        from pycombat import Combat
    else:
        print("Batch effect correction: None")
    
    # Set random seeds
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
        y = np.array(labels)
    
    if batch is not None:
        if isinstance(batch, pd.Series):
            batch_array = batch.values
        else:
            batch_array = np.array(batch)
    else:
        batch_array = None
    
    # Get all classes and encode labels
    classes = np.unique(y)
    n_classes = len(classes)
    print(f"Number of classes: {n_classes}")
    print("Class distribution:")
    for label in classes:
        count = np.sum(y == label)
        print(f"  {label}: {count}")
    
    # Encode labels to integers
    label_encoder = LabelEncoder()
    y_encoded = label_encoder.fit_transform(y)
    
    # Calculate class weights
    class_counts = np.bincount(y_encoded)
    total_samples = len(y_encoded)
    class_weights = total_samples / (n_classes * class_counts)
    class_weights_tensor = torch.FloatTensor(class_weights)
    print(f"Class weights: {dict(zip(classes, class_weights))}")
    
    if batch_array is not None:
        print("Batch distribution:")
        unique_batches, batch_counts = np.unique(
            batch_array,
            return_counts=True
        )
        for batch_name, count in zip(unique_batches, batch_counts):
            print(f"  {batch_name}: {count}")
    
    # Set device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # Initialize result storage
    n_total_folds = n_repeats * n_folds
    
    # Store metrics for each class
    results = {cls: {
        'auc': np.zeros(n_total_folds),
        'aupr': np.zeros(n_total_folds),
        'accuracy': np.zeros(n_total_folds),
        'recall': np.zeros(n_total_folds),
        'f1': np.zeros(n_total_folds)
    } for cls in classes}
    
    fold_counter = 0
    
    for rep in range(n_repeats):
        print(f"\n=== Repeat {rep + 1} ===")
        
        # Create stratified k-fold
        skf = StratifiedKFold(
            n_splits=n_folds,
            shuffle=True,
            random_state=random_state + rep
        )
        
        for fold, (train_idx, test_idx) in enumerate(skf.split(X, y_encoded)):
            print(f"  Fold {fold + 1}")
            
            # Split data
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]
            
            # Apply batch effect correction if requested
            if apply_combat and batch_array is not None:
                train_batch = batch_array[train_idx]
                test_batch = batch_array[test_idx]
                
                unique_train_batches = np.unique(train_batch)
                
                if len(unique_train_batches) >= 2:
                    unique_test_batches = np.unique(test_batch)
                    test_in_train = all(
                        b in unique_train_batches for b in unique_test_batches
                    )
                    
                    if test_in_train:
                        print(f"    Applying Combat correction "
                              f"({len(unique_train_batches)} batches)...")
                        
                        try:
                            combat = Combat()
                            X_train_corrected = combat.fit_transform(
                                Y=X_train,
                                b=train_batch
                            )
                            X_test_corrected = combat.transform(
                                Y=X_test,
                                b=test_batch
                            )
                            
                            if (np.isnan(X_train_corrected).any() or
                                    np.isinf(X_train_corrected).any()):
                                print("    Warning: Combat produced NaN/Inf, "
                                      "using original data")
                            else:
                                X_train = X_train_corrected
                                X_test = X_test_corrected
                                print("    Combat correction successful")
                        
                        except Exception as e:
                            print(f"    Error in ComBat: {str(e)}, "
                                  f"using original data")
                    else:
                        print("    Skipping ComBat: test batches not in "
                              "training set")
                else:
                    print(f"    Skipping ComBat: only "
                          f"{len(unique_train_batches)} batch in training set")
            
            print(f"    Training samples: {len(train_idx)}, "
                  f"Test samples: {len(test_idx)}")
            
            # Standardize features
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)
            
            # Convert to tensors
            X_train_tensor = torch.FloatTensor(X_train_scaled).to(device)
            y_train_tensor = torch.LongTensor(y_train).to(device)
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
            model = MultiClassNN(input_size, n_classes).to(device)
            
            # Loss function with class weights
            criterion = nn.CrossEntropyLoss(
                weight=class_weights_tensor.to(device)
            )
            optimizer = optim.Adam(model.parameters(), lr=lr)
            
            # Training
            model.train()
            for epoch in range(n_epochs):
                epoch_loss = 0.0
                for batch_X, batch_y in train_loader:
                    # Forward pass
                    outputs = model(batch_X)
                    loss = criterion(outputs, batch_y)
                    
                    # Backward pass
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                    
                    epoch_loss += loss.item()
            
            # Evaluation
            model.eval()
            with torch.no_grad():
                outputs = model(X_test_tensor)
                # Apply softmax to get probabilities
                y_pred_proba = torch.softmax(outputs, dim=1).cpu().numpy()
                y_pred = outputs.argmax(dim=1).cpu().numpy()
            
            # Overall accuracy
            acc_val = accuracy_score(y_test, y_pred)
            
            # Calculate per-class metrics using One-vs-Rest
            for cls in classes:
                # Get class index
                cls_idx = np.where(label_encoder.classes_ == cls)[0][0]
                
                # Binary labels for this class (OvR)
                y_test_binary = (y_test == cls_idx).astype(int)
                y_pred_binary = (y_pred == cls_idx).astype(int)
                
                # Get probability for this class
                y_pred_proba_cls = y_pred_proba[:, cls_idx]
                
                # AUC (One-vs-Rest)
                try:
                    auc_val = roc_auc_score(y_test_binary, y_pred_proba_cls)
                except Exception:
                    auc_val = np.nan
                
                # AUPR (One-vs-Rest)
                try:
                    aupr_val = average_precision_score(
                        y_test_binary,
                        y_pred_proba_cls
                    )
                except Exception:
                    aupr_val = np.nan
                
                # Recall and F1 for this class
                recall_val = recall_score(
                    y_test_binary,
                    y_pred_binary,
                    zero_division=0
                )
                f1_val = f1_score(
                    y_test_binary,
                    y_pred_binary,
                    zero_division=0
                )
                
                # Store results
                results[cls]['auc'][fold_counter] = auc_val
                results[cls]['aupr'][fold_counter] = aupr_val
                results[cls]['accuracy'][fold_counter] = acc_val
                results[cls]['recall'][fold_counter] = recall_val
                results[cls]['f1'][fold_counter] = f1_val
            
            fold_counter += 1
            
            print(f"    Overall Acc: {acc_val:.3f}")
            for cls in classes:
                print(f"    {cls} - "
                      f"AUC: {results[cls]['auc'][fold_counter-1]:.3f}, "
                      f"AUPR: {results[cls]['aupr'][fold_counter-1]:.3f}, "
                      f"Recall: {results[cls]['recall'][fold_counter-1]:.3f}, "
                      f"F1: {results[cls]['f1'][fold_counter-1]:.3f}")
    
    # Create metrics dataframe
    repeats_vec = np.repeat(np.arange(1, n_repeats + 1), n_folds)
    folds_vec = np.tile(np.arange(1, n_folds + 1), n_repeats)
    
    metrics_data = []
    for cls in classes:
        for metric in ['auc', 'aupr', 'accuracy', 'recall', 'f1']:
            for i in range(n_total_folds):
                metrics_data.append({
                    'disease': cls,
                    'metrics': metric,
                    'value': results[cls][metric][i],
                    'repeats': repeats_vec[i],
                    'folds': folds_vec[i]
                })
    
    metrics_df = pd.DataFrame(metrics_data)
    
    print(f"\nCross-validation completed!")
    print(f"Metrics dataframe shape: {metrics_df.shape}")
    print("\nOverall summary (per class):")
    summary = metrics_df.groupby(['disease', 'metrics'])['value'].agg(
        ['mean', 'std']
    )
    print(summary)
    
    return {
        'metrics_df': metrics_df
    }


def load_and_preprocess_data(otu_file, metadata_file, group_column='group',
                              sample_column='Sample', batch_column=None,
                              target_groups=None, abundance_threshold=0.001,
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
    batch_column : str or None
        Column name in metadata containing batch information
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
        - batch: Batch labels (if batch_column provided)
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
        raise ValueError(
            "No common samples found between OTU data and metadata"
        )
    
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
    
    # Get unique groups
    unique_groups = metadata[group_column].unique().tolist()
    
    # Create labels
    labels = metadata[group_column]
    
    print(f"\nClass distribution:")
    print(labels.value_counts())
    
    result = {
        'full_data': otu_filtered,
        'labels': labels,
        'metadata': metadata
    }
    
    # Get batch information if specified
    if batch_column is not None and batch_column in metadata.columns:
        batch = metadata[batch_column]
        result['batch'] = batch
        print(f"\nBatch distribution:")
        print(batch.value_counts())
    else:
        result['batch'] = None
    
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
    """Main function to run multi-disease prediction analysis."""
    parser = argparse.ArgumentParser(
        description='Multi-Disease Prediction using Multiple ML Methods'
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
    parser.add_argument(
        '--batch_column',
        type=str,
        default=None,
        help='Column name in metadata containing batch information '
             '(default: None)'
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
    
    # Batch correction
    parser.add_argument(
        '--apply_combat',
        action='store_true',
        help='Apply batch effect correction using pyComBat'
    )
    
    # Model selection
    parser.add_argument(
        '--methods',
        type=str,
        nargs='+',
        choices=['rf', 'nn', 'all'],
        default=['all'],
        help='Methods to run: rf, nn, or all (default: all)'
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
    parser.add_argument(
        '--animal',
        type=str,
        default=None,
        help='Animal type label to add to results (optional)'
    )
    
    args = parser.parse_args()
    
    # Determine which methods to run
    if 'all' in args.methods:
        methods_to_run = ['rf', 'nn']
    else:
        methods_to_run = args.methods
    
    print("="*70)
    print("Multi-Disease Prediction Analysis")
    print("="*70)
    print(f"Methods to run: {', '.join(methods_to_run)}")
    print(f"Cross-validation: {args.n_repeats} repeats x {args.n_folds} folds")
    print(f"Batch correction: {args.apply_combat}")
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
        batch_column=args.batch_column,
        target_groups=args.target_groups,
        abundance_threshold=args.abundance_threshold,
        signature_set=args.signature_set
    )
    
    full_data = data_dict['full_data']
    labels = data_dict['labels']
    batch = data_dict.get('batch', None)
    signature_data = data_dict.get('signature_data', None)
    
    results = {}
    
    # Run Random Forest
    if 'rf' in methods_to_run:
        print("\n" + "="*70)
        print("RANDOM FOREST - FULL DATASET")
        print("="*70)
        
        rf_full_results = rf_cv_multiclass(
            full_data,
            labels,
            batch=batch,
            apply_combat=args.apply_combat,
            n_trees=args.n_trees,
            n_folds=args.n_folds,
            n_repeats=args.n_repeats,
            random_state=args.random_state
        )
        results['rf_full'] = rf_full_results
        
        # Add metadata columns
        if args.animal is not None:
            rf_full_results['metrics_df']['animal'] = args.animal
        rf_full_results['metrics_df']['model'] = 'rf'
        rf_full_results['metrics_df']['feature'] = 'full'
        
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
            
            rf_sig_results = rf_cv_multiclass(
                signature_data,
                labels,
                batch=batch,
                apply_combat=args.apply_combat,
                n_trees=args.n_trees,
                n_folds=args.n_folds,
                n_repeats=args.n_repeats,
                random_state=args.random_state
            )
            results['rf_signature'] = rf_sig_results
            
            # Add metadata columns
            if args.animal is not None:
                rf_sig_results['metrics_df']['animal'] = args.animal
            rf_sig_results['metrics_df']['model'] = 'rf'
            rf_sig_results['metrics_df']['feature'] = 'signature'
            
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
    
    # Run Neural Network
    if 'nn' in methods_to_run:
        print("\n" + "="*70)
        print("NEURAL NETWORK - FULL DATASET")
        print("="*70)
        
        nn_full_results = nn_cv_multiclass(
            full_data,
            labels,
            batch=batch,
            apply_combat=args.apply_combat,
            n_epochs=args.n_epochs,
            batch_size=args.batch_size,
            lr=args.learning_rate,
            n_folds=args.n_folds,
            n_repeats=args.n_repeats,
            random_state=args.random_state
        )
        results['nn_full'] = nn_full_results
        
        # Add metadata columns
        if args.animal is not None:
            nn_full_results['metrics_df']['animal'] = args.animal
        nn_full_results['metrics_df']['model'] = 'network'
        nn_full_results['metrics_df']['feature'] = 'full'
        
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
            
            nn_sig_results = nn_cv_multiclass(
                signature_data,
                labels,
                batch=batch,
                apply_combat=args.apply_combat,
                n_epochs=args.n_epochs,
                batch_size=args.batch_size,
                lr=args.learning_rate,
                n_folds=args.n_folds,
                n_repeats=args.n_repeats,
                random_state=args.random_state
            )
            results['nn_signature'] = nn_sig_results
            
            # Add metadata columns
            if args.animal is not None:
                nn_sig_results['metrics_df']['animal'] = args.animal
            nn_sig_results['metrics_df']['model'] = 'network'
            nn_sig_results['metrics_df']['feature'] = 'signature'
            
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
        summary = result['metrics_df'].groupby(
            ['disease', 'metrics']
        )['value'].agg(['mean', 'std'])
        print(summary)
    
    print("\n" + "="*70)
    print("Analysis completed successfully!")
    print("="*70)
    
    return results


if __name__ == "__main__":
    main()

