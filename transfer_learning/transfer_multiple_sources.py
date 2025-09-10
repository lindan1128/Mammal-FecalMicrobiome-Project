import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, accuracy_score, recall_score, f1_score, confusion_matrix
from sklearn.utils import resample
from scipy import stats
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import random
import os


def set_seed(seed=42):
    """Set all random seeds for reproducibility"""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def extract_genus(taxon_names):
    """Extract genus names from full taxonomic strings"""
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


class SimpleNN(nn.Module):
    """Encoder Neural Network with freezing capabilities"""
    
    def __init__(self, input_size, encoding_dim=8):
        super(SimpleNN, self).__init__()
        # Encoder layers (dimensionality reduction)
        self.encoder1 = nn.Linear(input_size, 16)
        self.encoder2 = nn.Linear(16, encoding_dim)  # encoding dimension
        
        # Decoder layers (reconstruction) - for pre-training
        self.decoder1 = nn.Linear(encoding_dim, 16)
        self.decoder2 = nn.Linear(16, input_size)
        
        # Classification layer
        self.classifier = nn.Linear(encoding_dim, 1)
        
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        
    def forward(self, x, return_encoded=False):
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
        """Freeze encoder layers"""
        for param in self.encoder1.parameters():
            param.requires_grad = False
        for param in self.encoder2.parameters():
            param.requires_grad = False
    
    def unfreeze_encoder(self):
        """Unfreeze encoder layers"""
        for param in self.encoder1.parameters():
            param.requires_grad = True
        for param in self.encoder2.parameters():
            param.requires_grad = True
    
    def unfreeze_last_encoder_layer(self):
        """Unfreeze only the last encoder layer"""
        for param in self.encoder1.parameters():
            param.requires_grad = False
        for param in self.encoder2.parameters():
            param.requires_grad = True
    
    def get_encoded_features(self, x):
        """Get encoded features"""
        with torch.no_grad():
            encoded = self.relu(self.encoder1(x))
            encoded = self.relu(self.encoder2(encoded))
        return encoded


def pretrain_encoder(model, pretrain_data, pretrain_labels, epochs=50, 
                    reconstruction_weight=0.3, classification_weight=0.7):
    """Pre-train encoder on additional dataset (autoencoder + supervised classification)"""
    model.train()
    
    # Define loss functions
    mse_criterion = nn.MSELoss()  # reconstruction loss
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
        classification, encoded, decoded = model(pretrain_data, return_encoded=True)
        
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


def train_linear_probe(model, train_data, train_labels, epochs=50, lr=0.001):
    """Step-A: Freeze encoder, train only classification head"""
    model.freeze_encoder()
    model.train()
    
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.classifier.parameters(), lr=lr)
    
    for epoch in range(epochs):
        optimizer.zero_grad()
        outputs = model(train_data)
        loss = criterion(outputs, train_labels)
        loss.backward()
        optimizer.step()
    
    return model


def fine_tune_last_layer(model, train_data, train_labels, epochs=30, lr=0.0005):
    """Step-B: Unfreeze last encoder layer for fine-tuning"""
    model.unfreeze_last_encoder_layer()
    model.train()
    
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


def fine_tune_all_layers(model, train_data, train_labels, val_data, val_labels, 
                        epochs=50, lr=0.0001, patience=10):
    """Step-C: Unfreeze all layers for fine-tuning with early stopping"""
    model.unfreeze_encoder()
    model.train()
    
    criterion = nn.BCELoss()
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
            val_loss = criterion(val_outputs, val_labels)
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
                    epochs=100, lr=0.001, patience=15):
    """End-to-end training without freezing any layers"""
    model.unfreeze_encoder()  # Ensure all layers are trainable
    model.train()
    
    criterion = nn.BCELoss()
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
            val_loss = criterion(val_outputs, val_labels)
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
    """Evaluate model performance"""
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


def pretrain_encoder_multi_dataset(model, pretrain_data_list, pretrain_labels_list, epochs=50, 
                                  reconstruction_weight=0.3, classification_weight=0.7, batch_size=32):
    """Pre-train encoder on multiple datasets with balanced batch sampling"""
    model.train()
    
    # Define loss functions
    mse_criterion = nn.MSELoss()  # reconstruction loss
    bce_criterion = nn.BCELoss()  # classification loss
    
    # Optimizer includes all relevant parameters
    optimizer = optim.Adam([
        {'params': model.encoder1.parameters()},
        {'params': model.encoder2.parameters()},
        {'params': model.decoder1.parameters()},
        {'params': model.decoder2.parameters()},
        {'params': model.classifier.parameters()}
    ], lr=0.001)
    
    # Calculate samples per dataset in each batch
    num_datasets = len(pretrain_data_list)
    samples_per_dataset = batch_size // num_datasets
    remaining_samples = batch_size % num_datasets
    
    print(f"        Batch size: {batch_size}, Samples per dataset: {samples_per_dataset}")
    
    for epoch in range(epochs):
        epoch_total_loss = 0
        epoch_reconstruction_loss = 0
        epoch_classification_loss = 0
        num_batches = 0
        
        # Calculate number of batches per epoch (use minimum dataset size)
        min_dataset_size = min([len(data) for data in pretrain_data_list])
        num_batches_per_epoch = min_dataset_size // samples_per_dataset
        
        for batch_idx in range(num_batches_per_epoch):
            optimizer.zero_grad()
            
            # Sample from each dataset
            batch_data_list = []
            batch_labels_list = []
            
            for i, (data, labels) in enumerate(zip(pretrain_data_list, pretrain_labels_list)):
                # Calculate current dataset's batch size
                current_batch_size = samples_per_dataset
                if i < remaining_samples:  # Distribute remaining samples
                    current_batch_size += 1
                
                # Random sampling
                start_idx = (batch_idx * samples_per_dataset) % len(data)
                end_idx = start_idx + current_batch_size
                
                if end_idx > len(data):
                    # If exceeds range, continue from beginning
                    batch_data = torch.cat([data[start_idx:], data[:end_idx - len(data)]])
                    batch_labels = torch.cat([labels[start_idx:], labels[:end_idx - len(labels)]])
                else:
                    batch_data = data[start_idx:end_idx]
                    batch_labels = labels[start_idx:end_idx]
                
                batch_data_list.append(batch_data)
                batch_labels_list.append(batch_labels)
            
            # Combine all datasets' batches
            combined_batch_data = torch.cat(batch_data_list, dim=0)
            combined_batch_labels = torch.cat(batch_labels_list, dim=0)
            
            # Forward pass
            classification, encoded, decoded = model(combined_batch_data, return_encoded=True)
            
            # Calculate reconstruction loss
            reconstruction_loss = mse_criterion(decoded, combined_batch_data)
            
            # Calculate classification loss
            classification_loss = bce_criterion(classification, combined_batch_labels)
            
            # Combined loss
            total_loss = (reconstruction_weight * reconstruction_loss + 
                         classification_weight * classification_loss)
            
            # Backward pass
            total_loss.backward()
            optimizer.step()
            
            # Accumulate losses
            epoch_total_loss += total_loss.item()
            epoch_reconstruction_loss += reconstruction_loss.item()
            epoch_classification_loss += classification_loss.item()
            num_batches += 1
        
        # Print average losses every 10 epochs
        if (epoch + 1) % 10 == 0:
            avg_total_loss = epoch_total_loss / num_batches
            avg_reconstruction_loss = epoch_reconstruction_loss / num_batches
            avg_classification_loss = epoch_classification_loss / num_batches
            
            print(f"        Epoch {epoch+1}: Avg Total Loss: {avg_total_loss:.4f}, "
                  f"Avg Reconstruction: {avg_reconstruction_loss:.4f}, "
                  f"Avg Classification: {avg_classification_loss:.4f}")
    
    return model


def nn_cv_pretrain(data, labels, pretrain_data_list=None, pretrain_labels_list=None, n_folds=5, n_repeats=5, strategy="transfer"):
    """Cross-validation function with different training strategies"""
    if pretrain_data_list is not None:
        if strategy == "transfer":
            print(f"Running Neural Network with transfer learning: {n_folds} folds, {n_repeats} repeats")
            print("Training strategy: Pre-train encoder → Linear probe → Fine-tune all")
            print(f"Pre-training: Combined reconstruction + classification loss on {len(pretrain_data_list)} datasets")
        elif strategy == "cross":
            print(f"Running Neural Network with cross-dataset evaluation: {n_folds} folds, {n_repeats} repeats")
            print("Training strategy: Pre-train encoder → Direct evaluation (no fine-tuning)")
            print(f"Pre-training: Combined reconstruction + classification loss on {len(pretrain_data_list)} datasets")
    else:
        print(f"Running Neural Network with end-to-end training: {n_folds} folds, {n_repeats} repeats")
        print("Training strategy: End-to-end training (no pre-training)")
    
    # Store all metrics for each fold
    all_aucs = []
    all_accuracies = []
    all_recalls = []
    all_f1s = []
    all_fprs = []
    
    for rep in range(n_repeats):
        print(f"  Repeat {rep+1} - Balancing dataset...")
        
        # Balance the dataset
        normal_indices = np.where(labels == 0)[0]
        sick_indices = np.where(labels == 1)[0]
        
        if len(normal_indices) > len(sick_indices):
            # Downsample normal group
            balanced_normal = resample(normal_indices, n_samples=len(sick_indices), random_state=rep)
            balanced_indices = np.concatenate([balanced_normal, sick_indices])
        else:
            # Downsample sick group
            balanced_sick = resample(sick_indices, n_samples=len(normal_indices), random_state=rep)
            balanced_indices = np.concatenate([normal_indices, balanced_sick])
        
        balanced_data = data.iloc[balanced_indices]
        balanced_labels = labels[balanced_indices]
        
        print(f"    Original - Normal: {len(normal_indices)}, Sick: {len(sick_indices)}")
        print(f"    Balanced - Normal: {np.sum(balanced_labels == 0)}, Sick: {np.sum(balanced_labels == 1)}")
        
        # Create folds
        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=rep)
        
        for fold, (train_idx, test_idx) in enumerate(skf.split(balanced_data, balanced_labels)):
            print(f"    Fold {fold+1} - Training samples: {len(train_idx)}, Test samples: {len(test_idx)}")
            
            # Split data
            train_data = balanced_data.iloc[train_idx]
            train_labels = balanced_labels[train_idx]
            test_data = balanced_data.iloc[test_idx]
            test_labels = balanced_labels[test_idx]
            
            # Further split training set into train and validation (for early stopping)
            train_train_idx, val_idx = train_test_split(
                range(len(train_data)), test_size=0.2, 
                stratify=train_labels, random_state=rep
            )
            
            train_train_data = train_data.iloc[train_train_idx]
            train_train_labels = train_labels[train_train_idx]
            val_data = train_data.iloc[val_idx]
            val_labels = train_labels[val_idx]
            
            # Standardize data
            scaler = StandardScaler()
            train_train_scaled = scaler.fit_transform(train_train_data)
            val_scaled = scaler.transform(val_data)
            test_scaled = scaler.transform(test_data)
            
            # Convert to PyTorch tensors
            X_train_train = torch.FloatTensor(train_train_scaled)
            y_train_train = torch.FloatTensor(train_train_labels).unsqueeze(1)
            X_val = torch.FloatTensor(val_scaled)
            y_val = torch.FloatTensor(val_labels).unsqueeze(1)
            X_test = torch.FloatTensor(test_scaled)
            y_test = torch.FloatTensor(test_labels).unsqueeze(1)
            
            # Initialize model
            model = SimpleNN(input_size=train_data.shape[1])
            
            if pretrain_data_list is not None and pretrain_labels_list is not None:
                # With pre-training data: use step-wise training strategy
                # Combine pre-training data and training data for standardization
                all_pretrain_data = pd.concat(pretrain_data_list, ignore_index=True)
                combined_train_data = pd.concat([train_train_data, all_pretrain_data], ignore_index=True)
                print(f"      Combined training data: {len(combined_train_data)} samples")
                
                # Re-fit scaler with combined data
                scaler = StandardScaler()
                combined_train_scaled = scaler.fit_transform(combined_train_data)
                
                # Re-transform all data
                train_train_scaled = scaler.transform(train_train_data)
                val_scaled = scaler.transform(val_data)
                test_scaled = scaler.transform(test_data)
                
                # Re-convert to PyTorch tensors
                X_train_train = torch.FloatTensor(train_train_scaled)
                y_train_train = torch.FloatTensor(train_train_labels).unsqueeze(1)
                X_val = torch.FloatTensor(val_scaled)
                y_val = torch.FloatTensor(val_labels).unsqueeze(1)
                X_test = torch.FloatTensor(test_scaled)
                y_test = torch.FloatTensor(test_labels).unsqueeze(1)
                
                # Prepare pre-training data
                pretrain_tensors = []
                pretrain_label_tensors = []
                
                for i, (pretrain_data, pretrain_labels) in enumerate(zip(pretrain_data_list, pretrain_labels_list)):
                    pretrain_scaled = scaler.transform(pretrain_data)
                    X_pretrain = torch.FloatTensor(pretrain_scaled)
                    y_pretrain = torch.FloatTensor(pretrain_labels).unsqueeze(1)
                    pretrain_tensors.append(X_pretrain)
                    pretrain_label_tensors.append(y_pretrain)
                    print(f"      Dataset {i+1}: {len(pretrain_data)} samples")
                
                # Step 1: Pre-training encoder (multi-dataset, balanced batch)
                print(f"      Step 1: Pre-training encoder on {len(pretrain_data_list)} datasets...")
                model = pretrain_encoder_multi_dataset(model, pretrain_tensors, pretrain_label_tensors, 
                                                     epochs=30, reconstruction_weight=0.5, classification_weight=0.5)
                print(f"      Encoder pre-training completed")
                
                if strategy == "transfer":
                    # Step 2: Linear probe training (freeze encoder)
                    print(f"      Step 2: Linear probe training (frozen encoder)...")
                    model = train_linear_probe(model, X_train_train, y_train_train, epochs=100, lr=0.001)
                                    
                    # Step 3: Fine-tune all layers
                    print(f"      Step 3: Fine-tuning all layers...")
                    model = fine_tune_all_layers(model, X_train_train, y_train_train, X_val, y_val, epochs=100, lr=0.001)
                    final_auc, final_acc, final_rec, final_f1, final_fpr = evaluate_model(model, X_test, y_test)
                elif strategy == "cross":
                    # Cross strategy: directly evaluate on target dataset without fine-tuning
                    print(f"      Cross strategy: Direct evaluation on target dataset (no fine-tuning)...")
                    final_auc, final_acc, final_rec, final_f1, final_fpr = evaluate_model(model, X_test, y_test)
                    
            else:
                # No pre-training data: direct end-to-end training
                print(f"      End-to-end training (no pre-training)...")
                model = train_end_to_end(
                    model, X_train_train, y_train_train, X_val, y_val,
                    epochs=100, lr=0.001, patience=15
                )
                final_auc, final_acc, final_rec, final_f1, final_fpr = evaluate_model(model, X_test, y_test)
                print(f"      End-to-end training completed")
            
            # Store metrics
            all_aucs.append(final_auc)
            all_accuracies.append(final_acc)
            all_recalls.append(final_rec)
            all_f1s.append(final_f1)
            all_fprs.append(final_fpr)
            
            print(f"      Final - AUC: {final_auc:.3f}, Accuracy: {final_acc:.3f}, Recall: {final_rec:.3f}, F1: {final_f1:.3f}, FPR: {final_fpr:.3f}")
        
        # Calculate mean metrics for this repeat
        repeat_auc = np.mean(all_aucs[-n_folds:])
        repeat_accuracy = np.mean(all_accuracies[-n_folds:])
        repeat_recall = np.mean(all_recalls[-n_folds:])
        repeat_f1 = np.mean(all_f1s[-n_folds:])
        repeat_fpr = np.mean(all_fprs[-n_folds:])
        
        print(f"  Repeat {rep+1} completed - Mean AUC: {repeat_auc:.3f}, Accuracy: {repeat_accuracy:.3f}, Recall: {repeat_recall:.3f}, F1: {repeat_f1:.3f}, FPR: {repeat_fpr:.3f}")
    
    return {
        'aucs': np.array(all_aucs),
        'accuracies': np.array(all_accuracies),
        'recalls': np.array(all_recalls),
        'f1s': np.array(all_f1s),
        'fprs': np.array(all_fprs)
    }


def load_and_process_data(data_type="raw"):
    """Load and process cow datasets based on data type"""
    base_path = "/Users/lindan/Dropbox/PhD/Projects/Animal_Microbiome/OTU/cow/"
    
    # Define file suffixes based on data type
    if data_type == "raw":
        suffix = ""
    elif data_type == "mmuphin":
        suffix = "_MMUPHin"
    elif data_type == "limma":
        suffix = "_limma"
    else:
        raise ValueError("data_type must be 'raw', 'mmuphin', or 'limma'")
    
    # Load datasets
    datasets = {}
    for i, (otu_file, meta_file) in enumerate([
        ("25741" + suffix + ".csv", "PRJEB25741.csv"),
        ("716761" + suffix + ".csv", "PRJNA716761.csv"),
        ("918869" + suffix + ".csv", "PRJNA918869.csv")
    ], 1):
        print(f"Loading dataset {i}...")
        
        # Load OTU data
        otu_data = pd.read_csv(base_path + otu_file, index_col=0)
        
        # Load metadata
        metadata = pd.read_csv(base_path + meta_file, index_col=0)
        metadata = metadata[metadata['group'].isin(['normal', 'sick'])]
        
        # Find common samples
        metadata_samples = metadata['Sample'].tolist()
        common_samples = list(set(metadata_samples) & set(otu_data.columns))
        metadata = metadata[metadata['Sample'].isin(common_samples)]
        otu_data = otu_data[metadata['Sample']]
        
        # Reset indices
        metadata = metadata.reset_index(drop=True)
        
        # Process OTU data
        otu_transposed = otu_data.T
        otu_transposed.columns = extract_genus(otu_transposed.columns)
        
        # Filter low abundance taxa (< 0.001 in all samples)
        otu_filtered = otu_transposed.loc[:, (otu_transposed >= 0.001).any(axis=0)]
        
        # Handle duplicate columns by taking the mean
        otu_filtered = otu_filtered.groupby(level=0, axis=1).mean()
        
        # Process labels
        metadata['group_factor'] = metadata['group'].map({'normal': 0, 'sick': 1})
        
        datasets[i] = {
            'otu_data': otu_filtered,
            'metadata': metadata,
            'labels': metadata['group_factor'].values
        }
        
        print(f"  Dataset {i} - Samples: {otu_filtered.shape[0]}, Taxa: {otu_filtered.shape[1]}")
        print(f"  Group distribution: {metadata['group_factor'].value_counts().to_dict()}")
    
    return datasets


def align_datasets(datasets):
    """Align all datasets to have the same taxa"""
    # Find union of all taxa
    all_taxa = set()
    for dataset in datasets.values():
        all_taxa.update(dataset['otu_data'].columns)
    
    print(f"Total unique taxa across all datasets: {len(all_taxa)}")
    
    # Align each dataset
    aligned_datasets = {}
    for i, dataset in datasets.items():
        aligned_otu = dataset['otu_data'].reindex(columns=all_taxa, fill_value=0)
        aligned_datasets[i] = {
            'otu_data': aligned_otu,
            'metadata': dataset['metadata'],
            'labels': dataset['labels']
        }
        print(f"Dataset {i} aligned - Taxa: {aligned_otu.shape[1]}")
    
    return aligned_datasets


def get_signature_data(datasets):
    """Extract signature taxa from datasets"""
    signature_set = ['Barnesiella', 'Christensenellaceae_R-7_group', 'Odoribacter', 
                     '[Eubacterium]_coprostanoligenes_group', 'Agathobacter', 
                     'Bacteroidales_RF16_group', 'Megasphaera', 'Ruminococcus', 
                     'Alistipes', 'Parabacteroides', 'Phascolarctobacterium', 
                     'Succinivibrio', 'Oscillibacter', 'Enterococcus', 
                     'Clostridium_sensu_stricto_1', 'Escherichia-Shigella', 
                     'Prevotella', 'Prevotellaceae_Ga6A1_group', 'Prevotella_9', 'Collinsella']
    
    signature_datasets = {}
    for i, dataset in datasets.items():
        signature_data = dataset['otu_data'].loc[:, dataset['otu_data'].columns.isin(signature_set)]
        signature_datasets[i] = {
            'otu_data': signature_data,
            'metadata': dataset['metadata'],
            'labels': dataset['labels']
        }
        print(f"Dataset {i} signature taxa: {signature_data.shape[1]} / {len(signature_set)}")
    
    return signature_datasets


def run_combine_transfer_learning_analysis(datasets, feature_type="full"):
    """Run combine transfer learning analysis for all dataset combinations"""
    results = []
    
    # Select feature data
    if feature_type == "full":
        feature_datasets = datasets
    else:  # signature
        feature_datasets = get_signature_data(datasets)
    
    # Run analysis for each target dataset (using other two as pre-training)
    for target_id in [1, 2, 3]:
        print(f"\n=== Target Dataset {target_id} ===")
        target_data = feature_datasets[target_id]['otu_data']
        target_labels = feature_datasets[target_id]['labels']
        
        # Get source datasets (the other two)
        source_ids = [i for i in [1, 2, 3] if i != target_id]
        source_data_list = [feature_datasets[i]['otu_data'] for i in source_ids]
        source_labels_list = [feature_datasets[i]['labels'] for i in source_ids]
        
        print(f"Using datasets {source_ids} for pre-training")
        
        # Single dataset (no pre-training)
        print(f"Single dataset analysis (no pre-training):")
        single_results = nn_cv_pretrain(target_data, target_labels)
        
        # Store single dataset results
        for rep in range(5):
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
                    'model': 'single',
                    'target_dataset': target_id,
                    'source_dataset': None
                })
        
        # Multi dataset (with pre-training from other two datasets)
        # Transfer learning strategy
        print(f"Combine transfer learning from datasets {source_ids} to dataset {target_id}:")
        transfer_results = nn_cv_pretrain(
            target_data, target_labels,
            pretrain_data_list=source_data_list,
            pretrain_labels_list=source_labels_list,
            strategy="transfer"
        )
        
        # Store transfer learning results
        for rep in range(5):
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
                    'model': 'transfer',
                    'target_dataset': target_id,
                    'source_dataset': f"{source_ids[0]}+{source_ids[1]}"
                })
        
        # Cross strategy (pre-train on combined source datasets, evaluate on target without fine-tuning)
        print(f"Combine cross-dataset evaluation from datasets {source_ids} to dataset {target_id}:")
        cross_results = nn_cv_pretrain(
            target_data, target_labels,
            pretrain_data_list=source_data_list,
            pretrain_labels_list=source_labels_list,
            strategy="cross"
        )
        
        # Store cross strategy results
        for rep in range(5):
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
                    'model': 'cross',
                    'target_dataset': target_id,
                    'source_dataset': f"{source_ids[0]}+{source_ids[1]}"
                })
    
    return results


def main():
    """Main function to run the complete analysis"""
    # Set random seed
    set_seed(42)
    
    # Define data types to analyze
    data_types = ["raw", "mmuphin", "limma"]
    
    all_results = []
    
    for data_type in data_types:
        print(f"\n{'='*60}")
        print(f"ANALYZING DATA TYPE: {data_type.upper()}")
        print(f"{'='*60}")
        
        # Load and process data
        datasets = load_and_process_data(data_type)
        
        # Align datasets
        aligned_datasets = align_datasets(datasets)
        
        # Run analysis for full feature set
        print(f"\n--- Full Feature Set Analysis ---")
        full_results = run_combine_transfer_learning_analysis(aligned_datasets, "full")
        # Add data_type to each result
        for result in full_results:
            result['data_type'] = data_type
        all_results.extend(full_results)
        
        # Run analysis for signature set
        print(f"\n--- Signature Set Analysis ---")
        signature_results = run_combine_transfer_learning_analysis(aligned_datasets, "signature")
        # Add data_type to each result
        for result in signature_results:
            result['data_type'] = data_type
        all_results.extend(signature_results)
    
    # Convert results to DataFrame and save
    results_df = pd.DataFrame(all_results)
    
    # Save results
    output_file = "/Users/lindan/Dropbox/PhD/Projects/Animal_Microbiome/transfer_combine_cow.csv"
    results_df.to_csv(output_file, index=False)
    
    print(f"\n{'='*60}")
    print(f"ANALYSIS COMPLETED")
    print(f"{'='*60}")
    print(f"Results saved to: {output_file}")
    print(f"Total results: {len(results_df)} rows")
    print(f"Columns: {list(results_df.columns)}")
    
    # Print summary statistics
    print(f"\nSummary by model type:")
    summary = results_df.groupby(['model', 'feature']).agg({
        'AUC': ['mean', 'std'],
        'Acc': ['mean', 'std'],
        'Recall': ['mean', 'std'],
        'F1': ['mean', 'std'],
        'FPR': ['mean', 'std']
    }).round(3)
    print(summary)


if __name__ == "__main__":
    main()
