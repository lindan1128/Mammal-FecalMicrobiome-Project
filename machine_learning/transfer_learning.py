#!/usr/bin/env python3
"""
Transfer Learning for Disease Prediction

This script performs transfer learning analysis using neural networks:
- Single source transfer learning
- Multi-source transfer learning
- Cross-dataset evaluation
- With/without batch effect correction using pyComBat

The script supports multiple training strategies and evaluates the benefit
of transfer learning compared to single-dataset training.
"""

import argparse
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score,
    accuracy_score,
    recall_score,
    f1_score,
    confusion_matrix
)
from sklearn.utils import resample
from scipy import stats
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import random
import os
from pycombat import Combat


def set_seed(seed=42):
    """
    Set all random seeds for reproducibility.
    
    Parameters:
    -----------
    seed : int
        Random seed value
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


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
    genera = []
    for name in taxon_names:
        parts = name.split(';')
        if len(parts) >= 6:
            genus = parts[5]
            if genus.startswith('g__'):
                genus = genus[3:]
            genera.append(genus)
        else:
            genera.append(name)
    return genera


def load_and_process_data(otu_files, metadata_files, dataset_names,
                           group_column='group', sample_column='Sample',
                           target_groups=None, abundance_threshold=0.001):
    """
    Load and process multiple datasets with batch information.
    
    Parameters:
    -----------
    otu_files : list of str
        Paths to OTU table CSV files
    metadata_files : list of str
        Paths to metadata CSV files
    dataset_names : list of str
        Names for each dataset (used as batch labels)
    group_column : str
        Column name in metadata containing group labels
    sample_column : str
        Column name in metadata containing sample IDs
    target_groups : list or None
        List of target groups to include (e.g., ['normal', 'sick'])
    abundance_threshold : float
        Minimum abundance threshold for filtering taxa
    
    Returns:
    --------
    dict
        Dictionary mapping dataset index to data dictionary containing:
        - otu_data: Filtered OTU DataFrame
        - metadata: Metadata DataFrame
        - labels: Numeric labels (0/1)
        - batch: Batch labels (dataset names)
        - dataset_name: Name of the dataset
    """
    datasets = {}
    
    if target_groups is None:
        target_groups = ['normal', 'sick']
    
    # Create label mapping
    label_mapping = {target_groups[0]: 0, target_groups[1]: 1}
    
    for i, (otu_file, meta_file, dataset_name) in enumerate(
        zip(otu_files, metadata_files, dataset_names), 1
    ):
        print(f"Loading dataset {i} ({dataset_name})...")
        
        # Load OTU data
        otu_data = pd.read_csv(otu_file, index_col=0)
        
        # Load metadata
        metadata = pd.read_csv(meta_file, index_col=0)
        metadata = metadata[metadata[group_column].isin(target_groups)]
        
        # Find common samples
        metadata_samples = metadata[sample_column].tolist()
        common_samples = list(set(metadata_samples) & set(otu_data.columns))
        metadata = metadata[metadata[sample_column].isin(common_samples)]
        otu_data = otu_data[metadata[sample_column]]
        
        # Reset indices
        metadata = metadata.reset_index(drop=True)
        
        # Process OTU data
        otu_transposed = otu_data.T
        otu_transposed.columns = extract_genus(otu_transposed.columns)
        
        # Filter low abundance taxa
        otu_filtered = otu_transposed.loc[
            :, (otu_transposed >= abundance_threshold).any(axis=0)
        ]
        
        # Handle duplicate columns by taking the mean
        otu_filtered = otu_filtered.groupby(level=0, axis=1).mean()
        
        # Process labels
        metadata['group_factor'] = metadata[group_column].map(label_mapping)
        
        # Add batch information (dataset name as batch)
        metadata['batch'] = dataset_name
        
        datasets[i] = {
            'otu_data': otu_filtered,
            'metadata': metadata,
            'labels': metadata['group_factor'].values,
            'batch': metadata['batch'].values,
            'dataset_name': dataset_name
        }
        
        print(f"  Dataset {i} - Samples: {otu_filtered.shape[0]}, "
              f"Taxa: {otu_filtered.shape[1]}")
        print(f"  Group distribution: "
              f"{metadata['group_factor'].value_counts().to_dict()}")
    
    return datasets


def align_datasets(datasets):
    """
    Align all datasets to have the same taxa while preserving batch info.
    
    Parameters:
    -----------
    datasets : dict
        Dictionary of datasets from load_and_process_data
    
    Returns:
    --------
    dict
        Dictionary of aligned datasets with same taxa across all datasets
    """
    # Find union of all taxa
    all_taxa = set()
    for dataset in datasets.values():
        all_taxa.update(dataset['otu_data'].columns)
    
    print(f"Total unique taxa across all datasets: {len(all_taxa)}")
    
    # Align each dataset
    aligned_datasets = {}
    for i, dataset in datasets.items():
        aligned_otu = dataset['otu_data'].reindex(
            columns=all_taxa,
            fill_value=0
        )
        aligned_datasets[i] = {
            'otu_data': aligned_otu,
            'metadata': dataset['metadata'],
            'labels': dataset['labels'],
            'batch': dataset['batch'],
            'dataset_name': dataset['dataset_name']
        }
        print(f"Dataset {i} aligned - Taxa: {aligned_otu.shape[1]}")
    
    return aligned_datasets


def get_signature_data(datasets, signature_set):
    """
    Extract signature taxa from datasets while preserving batch information.
    
    Parameters:
    -----------
    datasets : dict
        Dictionary of datasets from align_datasets
    signature_set : list
        List of signature taxa names to extract
    
    Returns:
    --------
    dict
        Dictionary of datasets containing only signature taxa
    """
    signature_datasets = {}
    for i, dataset in datasets.items():
        signature_data = dataset['otu_data'].loc[
            :, dataset['otu_data'].columns.isin(signature_set)
        ]
        signature_datasets[i] = {
            'otu_data': signature_data,
            'metadata': dataset['metadata'],
            'labels': dataset['labels'],
            'batch': dataset['batch'],
            'dataset_name': dataset['dataset_name']
        }
        print(f"Dataset {i} signature taxa: {signature_data.shape[1]} / "
              f"{len(signature_set)}")
    
    return signature_datasets


class SimpleNN(nn.Module):
    """
    Encoder Neural Network with freezing capabilities for transfer learning.
    
    This network consists of:
    - Encoder layers (for feature learning and dimensionality reduction)
    - Decoder layers (for reconstruction during pre-training)
    - Classification layer (for final prediction)
    """
    
    def __init__(self, input_size, encoding_dim=8):
        """
        Initialize the neural network.
        
        Parameters:
        -----------
        input_size : int
            Number of input features
        encoding_dim : int
            Dimension of encoded representation
        """
        super(SimpleNN, self).__init__()
        # Encoder layers (dimensionality reduction)
        self.encoder1 = nn.Linear(input_size, 16)
        self.encoder2 = nn.Linear(16, encoding_dim)
        
        # Decoder layers (reconstruction) - for pre-training
        self.decoder1 = nn.Linear(encoding_dim, 16)
        self.decoder2 = nn.Linear(16, input_size)
        
        # Classification layer
        self.classifier = nn.Linear(encoding_dim, 1)
        
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        
    def forward(self, x, return_encoded=False):
        """
        Forward pass through the network.
        
        Parameters:
        -----------
        x : torch.Tensor
            Input data
        return_encoded : bool
            If True, return encoded representation and reconstruction
        
        Returns:
        --------
        torch.Tensor or tuple
            Classification output, or tuple of (classification, encoded,
            decoded) if return_encoded=True
        """
        # Encoder part
        encoded = self.relu(self.encoder1(x))
        encoded = self.relu(self.encoder2(encoded))
        
        # Decoder part (reconstruction) - only for pre-training
        decoded = self.relu(self.decoder1(encoded))
        decoded = self.sigmoid(self.decoder2(decoded))
        
        # Classification
        classification = self.sigmoid(self.classifier(encoded))
        
        if return_encoded:
            return classification, encoded, decoded
        return classification
    
    def freeze_encoder(self):
        """Freeze encoder layers to prevent weight updates."""
        for param in self.encoder1.parameters():
            param.requires_grad = False
        for param in self.encoder2.parameters():
            param.requires_grad = False
    
    def unfreeze_encoder(self):
        """Unfreeze encoder layers to allow weight updates."""
        for param in self.encoder1.parameters():
            param.requires_grad = True
        for param in self.encoder2.parameters():
            param.requires_grad = True
    
    def unfreeze_last_encoder_layer(self):
        """Unfreeze only the last encoder layer for gradual fine-tuning."""
        for param in self.encoder1.parameters():
            param.requires_grad = False
        for param in self.encoder2.parameters():
            param.requires_grad = True
    
    def get_encoded_features(self, x):
        """
        Get encoded features without gradient computation.
        
        Parameters:
        -----------
        x : torch.Tensor
            Input data
        
        Returns:
        --------
        torch.Tensor
            Encoded features
        """
        with torch.no_grad():
            encoded = self.relu(self.encoder1(x))
            encoded = self.relu(self.encoder2(encoded))
        return encoded


def pretrain_encoder(model, pretrain_data, pretrain_labels,
                      pretrain_weights=None, epochs=50,
                      reconstruction_weight=0.3, classification_weight=0.7):
    """
    Pre-train encoder using combined autoencoder and supervised classification.
    
    This function trains the encoder on source dataset(s) using a combination
    of reconstruction loss (autoencoder) and classification loss (supervised).
    
    Parameters:
    -----------
    model : SimpleNN
        Neural network model
    pretrain_data : torch.Tensor
        Pre-training feature data
    pretrain_labels : torch.Tensor
        Pre-training labels
    pretrain_weights : torch.Tensor or None
        Sample weights for weighted loss
    epochs : int
        Number of training epochs
    reconstruction_weight : float
        Weight for reconstruction loss in combined objective
    classification_weight : float
        Weight for classification loss in combined objective
    
    Returns:
    --------
    SimpleNN
        Pre-trained model
    """
    model.train()
    
    # Define loss functions
    mse_criterion = nn.MSELoss()  # reconstruction loss
    
    # Use weighted BCE loss if weights are provided
    if pretrain_weights is not None:
        bce_criterion = nn.BCELoss(weight=pretrain_weights)
    else:
        bce_criterion = nn.BCELoss()  # classification loss
    
    # Optimizer includes all relevant parameters
    optimizer = optim.Adam([
        {'params': model.encoder1.parameters()},
        {'params': model.encoder2.parameters()},
        {'params': model.decoder1.parameters()},
        {'params': model.decoder2.parameters()},
        {'params': model.classifier.parameters()}
    ], lr=0.001)
    
    for epoch in range(epochs):
        optimizer.zero_grad()
        
        # Forward pass
        classification, encoded, decoded = model(
            pretrain_data,
            return_encoded=True
        )
        
        # Calculate reconstruction loss
        reconstruction_loss = mse_criterion(decoded, pretrain_data)
        
        # Calculate classification loss
        classification_loss = bce_criterion(classification, pretrain_labels)
        
        # Combined loss
        total_loss = (reconstruction_weight * reconstruction_loss +
                      classification_weight * classification_loss)
        
        # Backward pass
        total_loss.backward()
        optimizer.step()
        
        # Print loss every 10 epochs
        if (epoch + 1) % 10 == 0:
            print(f"        Epoch {epoch+1}: Total Loss: {total_loss:.4f}, "
                  f"Reconstruction: {reconstruction_loss:.4f}, "
                  f"Classification: {classification_loss:.4f}")
    
    return model


def train_linear_probe(model, train_data, train_labels, sample_weights=None,
                        epochs=50, lr=0.001):
    """
    Step-A: Freeze encoder, train only classification head (linear probe).
    
    Parameters:
    -----------
    model : SimpleNN
        Neural network model with pre-trained encoder
    train_data : torch.Tensor
        Training feature data
    train_labels : torch.Tensor
        Training labels
    sample_weights : torch.Tensor or None
        Sample weights for weighted loss
    epochs : int
        Number of training epochs
    lr : float
        Learning rate
    
    Returns:
    --------
    SimpleNN
        Model with trained classification head
    """
    model.freeze_encoder()
    model.train()
    
    # Use weighted loss if weights are provided
    if sample_weights is not None:
        criterion = nn.BCELoss(weight=sample_weights)
    else:
        criterion = nn.BCELoss()
    optimizer = optim.Adam(model.classifier.parameters(), lr=lr)
    
    for epoch in range(epochs):
        optimizer.zero_grad()
        outputs = model(train_data)
        loss = criterion(outputs, train_labels)
        loss.backward()
        optimizer.step()
    
    return model


def fine_tune_last_layer(model, train_data, train_labels, sample_weights=None,
                          epochs=30, lr=0.001):
    """
    Step-B: Unfreeze last encoder layer for gradual fine-tuning.
    
    Parameters:
    -----------
    model : SimpleNN
        Neural network model
    train_data : torch.Tensor
        Training feature data
    train_labels : torch.Tensor
        Training labels
    sample_weights : torch.Tensor or None
        Sample weights for weighted loss
    epochs : int
        Number of training epochs
    lr : float
        Learning rate
    
    Returns:
    --------
    SimpleNN
        Model with fine-tuned last encoder layer
    """
    model.unfreeze_last_encoder_layer()
    model.train()
    
    # Use weighted loss if weights are provided
    if sample_weights is not None:
        criterion = nn.BCELoss(weight=sample_weights)
    else:
        criterion = nn.BCELoss()
    optimizer = optim.Adam([
        {'params': model.encoder2.parameters()},
        {'params': model.classifier.parameters()}
    ], lr=lr)
    
    for epoch in range(epochs):
        optimizer.zero_grad()
        outputs = model(train_data)
        loss = criterion(outputs, train_labels)
        loss.backward()
        optimizer.step()
    
    return model


def fine_tune_all_layers(model, train_data, train_labels, val_data,
                          val_labels, sample_weights=None, val_weights=None,
                          epochs=50, lr=0.0001, patience=10):
    """
    Step-C: Unfreeze all layers for full fine-tuning with early stopping.
    
    Parameters:
    -----------
    model : SimpleNN
        Neural network model
    train_data : torch.Tensor
        Training feature data
    train_labels : torch.Tensor
        Training labels
    val_data : torch.Tensor
        Validation feature data
    val_labels : torch.Tensor
        Validation labels
    sample_weights : torch.Tensor or None
        Sample weights for training loss
    val_weights : torch.Tensor or None
        Sample weights for validation loss
    epochs : int
        Maximum number of training epochs
    lr : float
        Learning rate
    patience : int
        Early stopping patience
    
    Returns:
    --------
    SimpleNN
        Fully fine-tuned model
    """
    model.unfreeze_encoder()
    model.train()
    
    # Use weighted loss if weights are provided
    if sample_weights is not None:
        criterion = nn.BCELoss(weight=sample_weights)
    else:
        criterion = nn.BCELoss()
    
    # Validation loss criterion
    if val_weights is not None:
        val_criterion = nn.BCELoss(weight=val_weights)
    else:
        val_criterion = nn.BCELoss()
    
    optimizer = optim.Adam(model.parameters(), lr=lr)
    
    best_val_loss = float('inf')
    patience_counter = 0
    
    for epoch in range(epochs):
        # Training
        optimizer.zero_grad()
        train_outputs = model(train_data)
        train_loss = criterion(train_outputs, train_labels)
        train_loss.backward()
        optimizer.step()
        
        # Validation
        model.eval()
        with torch.no_grad():
            val_outputs = model(val_data)
            val_loss = val_criterion(val_outputs, val_labels)
        model.train()
        
        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch+1}")
                break
    
    return model


def train_end_to_end(model, train_data, train_labels, val_data, val_labels,
                      sample_weights=None, val_weights=None, epochs=100,
                      lr=0.001, patience=15):
    """
    End-to-end training without freezing any layers (baseline approach).
    
    Parameters:
    -----------
    model : SimpleNN
        Neural network model
    train_data : torch.Tensor
        Training feature data
    train_labels : torch.Tensor
        Training labels
    val_data : torch.Tensor
        Validation feature data
    val_labels : torch.Tensor
        Validation labels
    sample_weights : torch.Tensor or None
        Sample weights for training loss
    val_weights : torch.Tensor or None
        Sample weights for validation loss
    epochs : int
        Maximum number of training epochs
    lr : float
        Learning rate
    patience : int
        Early stopping patience
    
    Returns:
    --------
    SimpleNN
        Trained model
    """
    model.unfreeze_encoder()  # Ensure all layers are trainable
    model.train()
    
    # Use weighted loss if weights are provided
    if sample_weights is not None:
        criterion = nn.BCELoss(weight=sample_weights)
    else:
        criterion = nn.BCELoss()
    
    # Validation loss criterion
    if val_weights is not None:
        val_criterion = nn.BCELoss(weight=val_weights)
    else:
        val_criterion = nn.BCELoss()
    
    optimizer = optim.Adam(model.parameters(), lr=lr)
    
    best_val_loss = float('inf')
    patience_counter = 0
    
    for epoch in range(epochs):
        # Training
        optimizer.zero_grad()
        train_outputs = model(train_data)
        train_loss = criterion(train_outputs, train_labels)
        train_loss.backward()
        optimizer.step()
        
        # Validation
        model.eval()
        with torch.no_grad():
            val_outputs = model(val_data)
            val_loss = val_criterion(val_outputs, val_labels)
        model.train()
        
        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch+1}")
                break
    
    return model


def evaluate_model(model, test_data, test_labels):
    """
    Evaluate model performance on test data.
    
    Parameters:
    -----------
    model : SimpleNN
        Trained neural network model
    test_data : torch.Tensor
        Test feature data
    test_labels : torch.Tensor or np.ndarray
        Test labels
    
    Returns:
    --------
    tuple
        (auc, accuracy, recall, f1, fpr)
    """
    model.eval()
    with torch.no_grad():
        test_probs = model(test_data).numpy().flatten()
        test_preds = (test_probs > 0.5).astype(int)
    
    # Calculate metrics
    auc = roc_auc_score(test_labels, test_probs)
    accuracy = accuracy_score(test_labels, test_preds)
    recall = recall_score(test_labels, test_preds)
    f1 = f1_score(test_labels, test_preds)
    
    # Calculate False Positive Rate
    tn, fp, fn, tp = confusion_matrix(test_labels, test_preds).ravel()
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
    
    return auc, accuracy, recall, f1, fpr


def calculate_sample_weights(labels):
    """
    Calculate sample weights for imbalanced datasets.
    
    Parameters:
    -----------
    labels : np.ndarray
        Array of labels (0 or 1)
    
    Returns:
    --------
    torch.Tensor
        Tensor of weights with same shape as labels
    """
    # Count samples for each class
    class_counts = np.bincount(labels.astype(int).flatten())
    n_samples = len(labels)
    
    # Calculate weight for each class (inverse of frequency)
    class_weights = n_samples / (len(class_counts) * class_counts)
    
    # Assign weight to each sample based on its class
    sample_weights = np.zeros(len(labels))
    for i, label in enumerate(labels.flatten()):
        sample_weights[i] = class_weights[int(label)]
    
    # Convert to torch tensor with same shape as labels
    return torch.FloatTensor(sample_weights).unsqueeze(1)


def apply_batch_correction(train_data, test_data, train_batch_labels,
                            test_batch_labels):
    """
    Apply batch effect correction using pyComBat with fit/transform pattern.
    
    This function uses the proper fit/transform pattern to avoid data leakage:
    1. Fit on training data to learn batch effect parameters
    2. Transform test data using parameters learned from training data
    
    Parameters:
    -----------
    train_data : pd.DataFrame
        Training samples (rows) and features (columns)
    test_data : pd.DataFrame
        Test samples (rows) and features (columns)
    train_batch_labels : pd.Series or np.ndarray
        Batch labels for training samples
    test_batch_labels : pd.Series or np.ndarray
        Batch labels for test samples
    
    Returns:
    --------
    tuple
        (train_corrected, test_corrected) - DataFrames with batch-corrected
        data
    """
    # Check if there are multiple batches
    unique_train_batches = pd.Series(train_batch_labels).nunique()
    
    # If only one batch exists in training data, skip correction
    if unique_train_batches < 2:
        print("      Only one batch detected in training data, "
              "skipping batch correction...")
        return train_data, test_data
    
    # Check if test batches are subset of train batches
    train_batch_set = set(pd.Series(train_batch_labels).unique())
    test_batch_set = set(pd.Series(test_batch_labels).unique())
    
    if not test_batch_set.issubset(train_batch_set):
        print(f"      Warning: Test batches {test_batch_set - train_batch_set} "
              f"not in training set")
        print("      Skipping batch correction to avoid errors...")
        return train_data, test_data
    
    # Use Combat class with fit/transform pattern
    combat = Combat()
    
    # Fit on training data
    print(f"      Fitting Combat on training data with "
          f"{unique_train_batches} batches...")
    train_corrected = combat.fit_transform(
        Y=train_data.values,
        b=train_batch_labels
    )
    
    # Transform test data
    print(f"      Transforming test data using training parameters...")
    test_corrected = combat.transform(
        Y=test_data.values,
        b=test_batch_labels
    )
    
    # Convert back to DataFrame
    train_corrected = pd.DataFrame(
        train_corrected,
        index=train_data.index,
        columns=train_data.columns
    )
    test_corrected = pd.DataFrame(
        test_corrected,
        index=test_data.index,
        columns=test_data.columns
    )
    
    return train_corrected, test_corrected


def nn_cv_pretrain(data, labels, batch_labels=None, pretrain_data=None,
                   pretrain_labels=None, pretrain_batch_labels=None,
                   n_folds=5, n_repeats=20, strategy="transfer",
                   apply_combat=False):
    """
    Cross-validation with transfer learning and different training strategies.
    
    Parameters:
    -----------
    data : pd.DataFrame
        OTU data for target dataset
    labels : np.ndarray
        Labels for target dataset
    batch_labels : pd.Series or np.ndarray or None
        Batch labels for target dataset
    pretrain_data : pd.DataFrame or None
        OTU data for source dataset(s)
    pretrain_labels : np.ndarray or None
        Labels for source dataset(s)
    pretrain_batch_labels : np.ndarray or None
        Batch labels for source dataset(s)
    n_folds : int
        Number of folds for cross-validation
    n_repeats : int
        Number of repeated cross-validation
    strategy : str
        Training strategy: "transfer" (pre-train + fine-tune) or
        "cross" (pre-train + direct evaluate)
    apply_combat : bool
        Whether to apply batch effect correction using pyComBat
    
    Returns:
    --------
    dict
        Dictionary containing arrays of:
        - aucs: AUC scores for each fold
        - accuracies: Accuracy scores
        - recalls: Recall scores
        - f1s: F1 scores
        - fprs: False positive rates
    """
    if pretrain_data is not None:
        if strategy == "transfer":
            print(f"Running Neural Network with transfer learning: "
                  f"{n_folds} folds, {n_repeats} repeats")
            print("Training strategy: Pre-train encoder → Linear probe → "
                  "Fine-tune all")
            print("Pre-training: Combined reconstruction + classification loss")
        elif strategy == "cross":
            print(f"Running Neural Network with cross-dataset evaluation: "
                  f"{n_folds} folds, {n_repeats} repeats")
            print("Training strategy: Pre-train encoder → "
                  "Direct evaluation (no fine-tuning)")
            print("Pre-training: Combined reconstruction + classification loss")
    else:
        print(f"Running Neural Network with end-to-end training: "
              f"{n_folds} folds, {n_repeats} repeats")
        print("Training strategy: End-to-end training (no pre-training)")
    
    if apply_combat:
        print("Batch effect correction: Using pyComBat within each fold")
    else:
        print("Batch effect correction: None")
    
    # Store all metrics for each fold
    all_aucs = []
    all_accuracies = []
    all_recalls = []
    all_f1s = []
    all_fprs = []
    
    for rep in range(n_repeats):
        print(f"  Repeat {rep+1} - Using full dataset with sample weights...")
        
        # No downsampling - use full dataset
        normal_count = np.sum(labels == 0)
        sick_count = np.sum(labels == 1)
        print(f"    Dataset - Normal: {normal_count}, Sick: {sick_count}")
        
        # Create folds
        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=rep)
        
        for fold, (train_idx, test_idx) in enumerate(skf.split(data, labels)):
            print(f"    Fold {fold+1} - Training samples: {len(train_idx)}, "
                  f"Test samples: {len(test_idx)}")
            
            # Split data
            train_data = data.iloc[train_idx]
            train_labels = labels[train_idx]
            test_data = data.iloc[test_idx]
            test_labels = labels[test_idx]
            
            # Get batch labels if provided
            if batch_labels is not None:
                train_batch = (batch_labels.iloc[train_idx]
                               if hasattr(batch_labels, 'iloc')
                               else batch_labels[train_idx])
                test_batch = (batch_labels.iloc[test_idx]
                              if hasattr(batch_labels, 'iloc')
                              else batch_labels[test_idx])
            else:
                train_batch = None
                test_batch = None
            
            # Further split training set into train and validation
            train_train_idx, val_idx = train_test_split(
                range(len(train_data)), test_size=0.2,
                stratify=train_labels, random_state=rep
            )
            
            train_train_data = train_data.iloc[train_train_idx]
            train_train_labels = train_labels[train_train_idx]
            val_data = train_data.iloc[val_idx]
            val_labels = train_labels[val_idx]
            
            # Calculate sample weights for imbalanced data
            train_train_weights = calculate_sample_weights(train_train_labels)
            val_weights = calculate_sample_weights(val_labels)
            
            # Initialize model
            model = SimpleNN(input_size=train_data.shape[1])
            
            if pretrain_data is not None and pretrain_labels is not None:
                # ========== Transfer Learning Path ==========
                # Combine pre-training data and training data
                combined_train_data = pd.concat(
                    [train_train_data, pretrain_data],
                    ignore_index=True
                )
                print(f"      Combined training data: "
                      f"{len(combined_train_data)} samples")
                
                # Apply batch effect correction if requested
                if apply_combat:
                    print(f"      Applying batch effect correction to "
                          f"combined data...")
                    
                    # Create batch labels for combined data
                    train_train_batch = (train_batch[train_train_idx]
                                         if train_batch is not None else None)
                    
                    if (train_train_batch is not None and
                            pretrain_batch_labels is not None):
                        # Combine batch labels
                        combined_batch = np.concatenate(
                            [train_train_batch, pretrain_batch_labels]
                        )
                        
                        # Check how many unique batches we have
                        unique_batches = pd.Series(combined_batch).nunique()
                        print(f"      Combined data has {unique_batches} "
                              f"unique batches")
                        
                        if unique_batches >= 2:
                            # Pre-check: remove low-variance features
                            print(f"      Pre-processing data for Combat...")
                            combined_data_array = combined_train_data.values
                            
                            # Calculate variance for each feature
                            feature_vars = np.var(combined_data_array, axis=0)
                            
                            # Identify features with sufficient variance
                            valid_features = feature_vars > 1e-10
                            n_valid = np.sum(valid_features)
                            n_removed = len(valid_features) - n_valid
                            
                            if n_removed > 0:
                                print(f"      Removing {n_removed} "
                                      f"low-variance features before Combat")
                                combined_data_filtered = (
                                    combined_data_array[:, valid_features]
                                )
                            else:
                                combined_data_filtered = combined_data_array
                            
                            # Apply Combat only if we have enough features
                            if n_valid > 0:
                                combat = Combat()
                                
                                print(f"      Fitting Combat on {n_valid} "
                                      f"features...")
                                try:
                                    combined_train_corrected_filtered = (
                                        combat.fit_transform(
                                            Y=combined_data_filtered,
                                            b=combined_batch
                                        )
                                    )
                                    
                                    # Check for NaN or Inf values
                                    n_nan = np.isnan(
                                        combined_train_corrected_filtered
                                    ).sum()
                                    n_inf = np.isinf(
                                        combined_train_corrected_filtered
                                    ).sum()
                                    
                                    if n_nan > 0 or n_inf > 0:
                                        print(f"      Warning: Combat produced "
                                              f"{n_nan} NaN and {n_inf} Inf "
                                              f"values")
                                        print(f"      Skipping Combat due to "
                                              f"numerical issues")
                                        combat_success = False
                                    else:
                                        print(f"      Combat succeeded without "
                                              f"numerical issues")
                                        combat_success = True
                                        
                                        # Reconstruct full array
                                        combined_train_corrected = np.zeros_like(
                                            combined_data_array
                                        )
                                        combined_train_corrected[
                                            :, valid_features
                                        ] = combined_train_corrected_filtered
                                        combined_train_corrected[
                                            :, ~valid_features
                                        ] = combined_data_array[
                                            :, ~valid_features
                                        ]
                                        
                                except Exception as e:
                                    print(f"      Error in Combat: {str(e)}")
                                    print(f"      Skipping Combat")
                                    combat_success = False
                            else:
                                print(f"      No valid features for Combat, "
                                      f"skipping")
                                combat_success = False
                            
                            if combat_success:
                                # Extract corrected portions
                                n_train_train = len(train_train_data)
                                train_train_data_corrected = (
                                    combined_train_corrected[:n_train_train]
                                )
                                pretrain_data_corrected = (
                                    combined_train_corrected[n_train_train:]
                                )
                                
                                # Convert back to DataFrame
                                train_train_data = pd.DataFrame(
                                    train_train_data_corrected,
                                    index=train_train_data.index,
                                    columns=train_train_data.columns
                                )
                                pretrain_data = pd.DataFrame(
                                    pretrain_data_corrected,
                                    index=pretrain_data.index,
                                    columns=pretrain_data.columns
                                )
                                
                                print(f"      Batch correction completed "
                                      f"successfully")
                            else:
                                print(f"      Using original data "
                                      f"(Combat skipped)")
                        else:
                            print(f"      Only one batch in combined data, "
                                  f"skipping batch correction")
                
                # ========== Standardization with Combined Data ==========
                print(f"      Standardizing data using combined dataset "
                      f"statistics...")
                combined_train_data = pd.concat(
                    [train_train_data, pretrain_data],
                    ignore_index=True
                )
                scaler = StandardScaler()
                combined_train_scaled = scaler.fit_transform(
                    combined_train_data
                )
                
                # Transform all datasets using the same scaler
                train_train_scaled = scaler.transform(train_train_data)
                val_scaled = scaler.transform(val_data)
                test_scaled = scaler.transform(test_data)
                pretrain_scaled = scaler.transform(pretrain_data)
                
                # Convert to PyTorch tensors
                X_train_train = torch.FloatTensor(train_train_scaled)
                y_train_train = torch.FloatTensor(
                    train_train_labels
                ).unsqueeze(1)
                X_val = torch.FloatTensor(val_scaled)
                y_val = torch.FloatTensor(val_labels).unsqueeze(1)
                X_test = torch.FloatTensor(test_scaled)
                y_test = torch.FloatTensor(test_labels).unsqueeze(1)
                X_pretrain = torch.FloatTensor(pretrain_scaled)
                y_pretrain = torch.FloatTensor(pretrain_labels).unsqueeze(1)
                
                # Calculate pretrain weights
                pretrain_weights = calculate_sample_weights(pretrain_labels)
                
                # ========== Step 1: Pre-training encoder ==========
                print(f"      Step 1: Pre-training encoder on "
                      f"{len(pretrain_data)} samples...")
                model = pretrain_encoder(
                    model, X_pretrain, y_pretrain, pretrain_weights,
                    epochs=100, reconstruction_weight=0.5,
                    classification_weight=0.5
                )
                print(f"      Encoder pre-training completed")
                
                if strategy == "transfer":
                    # Step 2: Linear probe training (freeze encoder)
                    print(f"      Step 2: Linear probe training "
                          f"(frozen encoder)...")
                    model = train_linear_probe(
                        model, X_train_train, y_train_train,
                        train_train_weights, epochs=100, lr=0.001
                    )
                    
                    # Step 3: Fine-tune all layers
                    print(f"      Step 3: Fine-tuning all layers...")
                    model = fine_tune_all_layers(
                        model, X_train_train, y_train_train, X_val, y_val,
                        train_train_weights, val_weights, epochs=100, lr=0.0001
                    )
                    final_auc, final_acc, final_rec, final_f1, final_fpr = (
                        evaluate_model(model, X_test, y_test)
                    )
                
                elif strategy == "cross":
                    # Cross strategy: evaluate without fine-tuning
                    print(f"      Cross strategy: Direct evaluation on target "
                          f"dataset (no fine-tuning)...")
                    final_auc, final_acc, final_rec, final_f1, final_fpr = (
                        evaluate_model(model, X_test, y_test)
                    )
                
            else:
                # ========== Direct Training Path (No Transfer Learning) ==========
                print(f"      Standardizing data using target dataset "
                      f"statistics...")
                scaler = StandardScaler()
                train_train_scaled = scaler.fit_transform(train_train_data)
                val_scaled = scaler.transform(val_data)
                test_scaled = scaler.transform(test_data)
                
                # Convert to PyTorch tensors
                X_train_train = torch.FloatTensor(train_train_scaled)
                y_train_train = torch.FloatTensor(
                    train_train_labels
                ).unsqueeze(1)
                X_val = torch.FloatTensor(val_scaled)
                y_val = torch.FloatTensor(val_labels).unsqueeze(1)
                X_test = torch.FloatTensor(test_scaled)
                y_test = torch.FloatTensor(test_labels).unsqueeze(1)
                
                # ========== End-to-end training ==========
                print(f"      End-to-end training (no pre-training)...")
                model = train_end_to_end(
                    model, X_train_train, y_train_train, X_val, y_val,
                    train_train_weights, val_weights, epochs=100, lr=0.001,
                    patience=15
                )
                final_auc, final_acc, final_rec, final_f1, final_fpr = (
                    evaluate_model(model, X_test, y_test)
                )
                print(f"      End-to-end training completed")
            
            # Store metrics
            all_aucs.append(final_auc)
            all_accuracies.append(final_acc)
            all_recalls.append(final_rec)
            all_f1s.append(final_f1)
            all_fprs.append(final_fpr)
            
            print(f"      Final - AUC: {final_auc:.3f}, "
                  f"Accuracy: {final_acc:.3f}, Recall: {final_rec:.3f}, "
                  f"F1: {final_f1:.3f}, FPR: {final_fpr:.3f}")
        
        # Calculate mean metrics for this repeat
        repeat_auc = np.mean(all_aucs[-n_folds:])
        repeat_accuracy = np.mean(all_accuracies[-n_folds:])
        repeat_recall = np.mean(all_recalls[-n_folds:])
        repeat_f1 = np.mean(all_f1s[-n_folds:])
        repeat_fpr = np.mean(all_fprs[-n_folds:])
        
        print(f"  Repeat {rep+1} completed - Mean AUC: {repeat_auc:.3f}, "
              f"Accuracy: {repeat_accuracy:.3f}, Recall: {repeat_recall:.3f}, "
              f"F1: {repeat_f1:.3f}, FPR: {repeat_fpr:.3f}")
    
    return {
        'aucs': np.array(all_aucs),
        'accuracies': np.array(all_accuracies),
        'recalls': np.array(all_recalls),
        'f1s': np.array(all_f1s),
        'fprs': np.array(all_fprs)
    }


def run_transfer_learning_analysis(datasets, feature_type="full",
                                     apply_combat=False, combat_suffix="",
                                     signature_set=None):
    """
    Run transfer learning analysis for all dataset combinations.
    
    Parameters:
    -----------
    datasets : dict
        Dictionary of datasets
    feature_type : str
        "full" or "signature"
    apply_combat : bool
        Whether to apply batch effect correction
    combat_suffix : str
        Suffix to add to model type if combat is applied
    signature_set : list or None
        List of signature taxa (required if feature_type="signature")
    
    Returns:
    --------
    list
        List of result dictionaries
    """
    results = []
    
    # Select feature data
    if feature_type == "full":
        feature_datasets = datasets
    else:  # signature
        if signature_set is None:
            raise ValueError(
                "signature_set must be provided when feature_type='signature'"
            )
        feature_datasets = get_signature_data(datasets, signature_set)
    
    # Run analysis for each target dataset
    for target_id in sorted(datasets.keys()):
        print(f"\n=== Target Dataset {target_id} ===")
        target_data = feature_datasets[target_id]['otu_data']
        target_labels = feature_datasets[target_id]['labels']
        target_batch = feature_datasets[target_id]['batch']
        
        # Single dataset (no pre-training)
        print(f"Single dataset analysis (no pre-training):")
        single_results = nn_cv_pretrain(
            target_data, target_labels, target_batch,
            apply_combat=apply_combat
        )
        
        # Store single dataset results
        for rep in range(20):
            for fold in range(5):
                idx = rep * 5 + fold
                results.append({
                    'AUC': single_results['aucs'][idx],
                    'Acc': single_results['accuracies'][idx],
                    'Recall': single_results['recalls'][idx],
                    'F1': single_results['f1s'][idx],
                    'FPR': single_results['fprs'][idx],
                    'repeats': rep + 1,
                    'fold': fold + 1,
                    'feature': feature_type,
                    'model': 'single' + combat_suffix,
                    'target_dataset': target_id,
                    'source_dataset': None,
                    'num_sources': 0,
                    'batch_correction': 'combat' if apply_combat else 'none'
                })
        
        # Multi dataset (with pre-training)
        # 1. Single source transfer learning
        for source_id in sorted(datasets.keys()):
            if source_id != target_id:
                source_data = feature_datasets[source_id]['otu_data']
                source_labels = feature_datasets[source_id]['labels']
                source_batch = feature_datasets[source_id]['batch']
                
                # Transfer learning strategy
                print(f"Transfer learning from dataset {source_id} to "
                      f"dataset {target_id}:")
                transfer_results = nn_cv_pretrain(
                    target_data, target_labels, target_batch,
                    pretrain_data=source_data,
                    pretrain_labels=source_labels,
                    pretrain_batch_labels=source_batch,
                    strategy="transfer",
                    apply_combat=apply_combat
                )
                
                # Store transfer learning results
                for rep in range(20):
                    for fold in range(5):
                        idx = rep * 5 + fold
                        results.append({
                            'AUC': transfer_results['aucs'][idx],
                            'Acc': transfer_results['accuracies'][idx],
                            'Recall': transfer_results['recalls'][idx],
                            'F1': transfer_results['f1s'][idx],
                            'FPR': transfer_results['fprs'][idx],
                            'repeats': rep + 1,
                            'fold': fold + 1,
                            'feature': feature_type,
                            'model': 'transfer' + combat_suffix,
                            'target_dataset': target_id,
                            'source_dataset': str(source_id),
                            'num_sources': 1,
                            'batch_correction': 'combat' if apply_combat else 'none'
                        })
                
                # Cross strategy
                print(f"Cross-dataset evaluation from dataset {source_id} to "
                      f"dataset {target_id}:")
                cross_results = nn_cv_pretrain(
                    target_data, target_labels, target_batch,
                    pretrain_data=source_data,
                    pretrain_labels=source_labels,
                    pretrain_batch_labels=source_batch,
                    strategy="cross",
                    apply_combat=apply_combat
                )
                
                # Store cross strategy results
                for rep in range(20):
                    for fold in range(5):
                        idx = rep * 5 + fold
                        results.append({
                            'AUC': cross_results['aucs'][idx],
                            'Acc': cross_results['accuracies'][idx],
                            'Recall': cross_results['recalls'][idx],
                            'F1': cross_results['f1s'][idx],
                            'FPR': cross_results['fprs'][idx],
                            'repeats': rep + 1,
                            'fold': fold + 1,
                            'feature': feature_type,
                            'model': 'cross' + combat_suffix,
                            'target_dataset': target_id,
                            'source_dataset': str(source_id),
                            'num_sources': 1,
                            'batch_correction': 'combat' if apply_combat else 'none'
                        })
        
        # 2. Multi-source transfer learning
        all_sources = sorted(datasets.keys())
        source_ids = [s for s in all_sources if s != target_id]
        
        if len(source_ids) >= 2:
            # Combine all sources
            print(f"\nMulti-source transfer learning from datasets "
                  f"{source_ids} to dataset {target_id}:")
            
            # Concatenate data from all sources
            multi_source_data_list = []
            multi_source_labels_list = []
            multi_source_batch_list = []
            
            for sid in source_ids:
                multi_source_data_list.append(
                    feature_datasets[sid]['otu_data']
                )
                multi_source_labels_list.append(
                    feature_datasets[sid]['labels']
                )
                multi_source_batch_list.append(feature_datasets[sid]['batch'])
            
            # Combine all sources
            multi_source_data = pd.concat(
                multi_source_data_list,
                ignore_index=True
            )
            multi_source_labels = np.concatenate(multi_source_labels_list)
            multi_source_batch = np.concatenate(multi_source_batch_list)
            
            print(f"  Combined source data: {len(multi_source_data)} samples "
                  f"from {len(source_ids)} datasets")
            
            # Transfer learning with multi-source
            multi_transfer_results = nn_cv_pretrain(
                target_data, target_labels, target_batch,
                pretrain_data=multi_source_data,
                pretrain_labels=multi_source_labels,
                pretrain_batch_labels=multi_source_batch,
                strategy="transfer",
                apply_combat=apply_combat
            )
            
            # Store multi-source transfer learning results
            source_str = ','.join(map(str, source_ids))
            for rep in range(20):
                for fold in range(5):
                    idx = rep * 5 + fold
                    results.append({
                        'AUC': multi_transfer_results['aucs'][idx],
                        'Acc': multi_transfer_results['accuracies'][idx],
                        'Recall': multi_transfer_results['recalls'][idx],
                        'F1': multi_transfer_results['f1s'][idx],
                        'FPR': multi_transfer_results['fprs'][idx],
                        'repeats': rep + 1,
                        'fold': fold + 1,
                        'feature': feature_type,
                        'model': 'transfer_multi' + combat_suffix,
                        'target_dataset': target_id,
                        'source_dataset': source_str,
                        'num_sources': len(source_ids),
                        'batch_correction': 'combat' if apply_combat else 'none'
                    })
            
            # Cross strategy with multi-source
            print(f"Cross-dataset evaluation from datasets {source_ids} to "
                  f"dataset {target_id}:")
            multi_cross_results = nn_cv_pretrain(
                target_data, target_labels, target_batch,
                pretrain_data=multi_source_data,
                pretrain_labels=multi_source_labels,
                pretrain_batch_labels=multi_source_batch,
                strategy="cross",
                apply_combat=apply_combat
            )
            
            # Store multi-source cross strategy results
            for rep in range(20):
                for fold in range(5):
                    idx = rep * 5 + fold
                    results.append({
                        'AUC': multi_cross_results['aucs'][idx],
                        'Acc': multi_cross_results['accuracies'][idx],
                        'Recall': multi_cross_results['recalls'][idx],
                        'F1': multi_cross_results['f1s'][idx],
                        'FPR': multi_cross_results['fprs'][idx],
                        'repeats': rep + 1,
                        'fold': fold + 1,
                        'feature': feature_type,
                        'model': 'cross_multi' + combat_suffix,
                        'target_dataset': target_id,
                        'source_dataset': source_str,
                        'num_sources': len(source_ids),
                        'batch_correction': 'combat' if apply_combat else 'none'
                    })
    
    return results


def main():
    """Main function to run transfer learning analysis."""
    parser = argparse.ArgumentParser(
        description='Transfer Learning for Disease Prediction'
    )
    
    # Input files
    parser.add_argument(
        '--otu_files',
        type=str,
        nargs='+',
        required=True,
        help='Paths to OTU table CSV files for each dataset'
    )
    parser.add_argument(
        '--metadata_files',
        type=str,
        nargs='+',
        required=True,
        help='Paths to metadata CSV files for each dataset'
    )
    parser.add_argument(
        '--dataset_names',
        type=str,
        nargs='+',
        required=True,
        help='Names for each dataset (used as batch labels)'
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
        default=['normal', 'sick'],
        help='Target groups to include in analysis (default: normal sick)'
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
        help='List of signature taxa to use (optional)'
    )
    
    # Analysis options
    parser.add_argument(
        '--feature_types',
        type=str,
        nargs='+',
        choices=['full', 'signature', 'both'],
        default=['both'],
        help='Feature types to analyze: full, signature, or both '
             '(default: both)'
    )
    parser.add_argument(
        '--batch_correction',
        type=str,
        nargs='+',
        choices=['none', 'combat', 'both'],
        default=['both'],
        help='Batch correction strategies: none, combat, or both '
             '(default: both)'
    )
    
    # Output options
    parser.add_argument(
        '--output_file',
        type=str,
        required=True,
        help='Path to output CSV file for results'
    )
    parser.add_argument(
        '--random_seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    
    args = parser.parse_args()
    
    # Validate input lengths
    if not (len(args.otu_files) == len(args.metadata_files) ==
            len(args.dataset_names)):
        raise ValueError(
            "Number of OTU files, metadata files, and dataset names must match"
        )
    
    # Determine feature types to run
    if 'both' in args.feature_types:
        feature_types = ['full', 'signature']
    else:
        feature_types = args.feature_types
    
    # Determine batch correction strategies
    if 'both' in args.batch_correction:
        batch_strategies = [(False, ""), (True, "_combat")]
    elif 'combat' in args.batch_correction:
        batch_strategies = [(True, "_combat")]
    else:
        batch_strategies = [(False, "")]
    
    print("="*70)
    print("Transfer Learning Analysis")
    print("="*70)
    print(f"Number of datasets: {len(args.dataset_names)}")
    print(f"Dataset names: {', '.join(args.dataset_names)}")
    print(f"Feature types: {', '.join(feature_types)}")
    print(f"Batch correction: "
          f"{', '.join([s[1] if s[1] else 'none' for s in batch_strategies])}")
    print("="*70)
    
    # Set random seed
    set_seed(args.random_seed)
    
    all_results = []
    
    # Load and process data
    print(f"\n{'='*70}")
    print(f"LOADING DATA")
    print(f"{'='*70}")
    datasets = load_and_process_data(
        args.otu_files,
        args.metadata_files,
        args.dataset_names,
        group_column=args.group_column,
        sample_column=args.sample_column,
        target_groups=args.target_groups,
        abundance_threshold=args.abundance_threshold
    )
    
    # Align datasets
    aligned_datasets = align_datasets(datasets)
    
    # Run analyses
    for apply_combat, suffix in batch_strategies:
        print(f"\n{'='*70}")
        if apply_combat:
            print(f"RUNNING ANALYSIS WITH PYCOMBAT BATCH CORRECTION")
        else:
            print(f"RUNNING ANALYSIS WITHOUT BATCH CORRECTION")
        print(f"{'='*70}")
        
        for feature_type in feature_types:
            if feature_type == 'signature' and args.signature_set is None:
                print(f"\nSkipping signature analysis: "
                      f"no signature set provided")
                continue
            
            print(f"\n--- {feature_type.capitalize()} Feature Set Analysis ---")
            results = run_transfer_learning_analysis(
                aligned_datasets,
                feature_type,
                apply_combat=apply_combat,
                combat_suffix=suffix,
                signature_set=args.signature_set
            )
            all_results.extend(results)
    
    # Convert results to DataFrame and save
    results_df = pd.DataFrame(all_results)
    
    # Save results
    results_df.to_csv(args.output_file, index=False)
    
    print(f"\n{'='*70}")
    print(f"ANALYSIS COMPLETED")
    print(f"{'='*70}")
    print(f"Results saved to: {args.output_file}")
    print(f"Total results: {len(results_df)} rows")
    print(f"Columns: {list(results_df.columns)}")
    
    # Print summary statistics
    print(f"\nSummary by model type and batch correction:")
    summary = results_df.groupby(
        ['model', 'feature', 'batch_correction']
    ).agg({
        'AUC': ['mean', 'std'],
        'Acc': ['mean', 'std'],
        'Recall': ['mean', 'std'],
        'F1': ['mean', 'std'],
        'FPR': ['mean', 'std']
    }).round(3)
    print(summary)


if __name__ == "__main__":
    main()
