import os
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score, recall_score
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# 1. Data loading
DATA_DIR = ''

profile_files = {
    'cow1': 'df_sub_core_cow1_corrected.csv',
    'cow2': 'df_sub_core_cow2_corrected.csv',
    'pig1': 'df_sub_core_pig1_corrected.csv',
    'pig2': 'df_sub_core_pig2_corrected.csv',
    'dog1': 'df_sub_core_dog1_corrected.csv',
    'edd_singh': 'df_sub_core_edd_singh_corrected.csv',
    'ibd_alm': 'df_sub_core_ibd_alm_corrected.csv',
    'ibd_engstrand': 'df_sub_core_ibd_engstrand_corrected.csv',
    'ibd_gevers': 'df_sub_core_ibd_gevers_corrected.csv',
    'ibd_huttenhower': 'df_sub_core_ibd_huttenhower_corrected.csv',
    'cdi_schubert': 'df_sub_core_cdi_schubert_corrected.csv',
    'cdi_vincent': 'df_sub_core_cdi_vincent_corrected.csv',
    'cdi_youngster': 'df_sub_core_cdi_youngster_corrected.csv'
}

label_files = {k: f'group_{k}.csv' for k in profile_files}

profiles = {k: pd.read_csv(os.path.join(DATA_DIR, v), index_col=0).T for k, v in profile_files.items()}
labels = {k: pd.read_csv(os.path.join(DATA_DIR, label_files[k]))['label'].values for k in profile_files}

# 2. Host grouping
host_map = {
    'cow': ['cow1', 'cow2'],
    'pig': ['pig1', 'pig2'],
    'dog': ['dog1'],
    'human': ['edd_singh', 'ibd_alm', 'ibd_engstrand', 'ibd_gevers', 'ibd_huttenhower', 'cdi_schubert', 'cdi_vincent', 'cdi_youngster']
}

# 3. DNN model definitions
class SimpleDNN(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 80),
            nn.ReLU(),
            nn.Linear(80, 40),
            nn.ReLU(),
            nn.Linear(40, 20),
            nn.ReLU(),
            nn.Linear(20, 1)
        )
    def forward(self, x):
        return self.net(x).squeeze(-1)

class AttentionMLP(nn.Module):
    def __init__(self, input_dim, hidden_dim=48, dropout_rate=0.3):
        super().__init__()
        self.hidden_dim = hidden_dim
        
        # Feature projection layer
        self.input_proj = nn.Linear(input_dim, hidden_dim)
        
        # Attention mechanism
        self.attention = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.Tanh(),
            nn.Linear(hidden_dim//2, 1)
        )
        
        # Feedforward network
        self.feed_forward = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.ReLU(),
            nn.Dropout(dropout_rate)
        )
        
        # Output layer
        self.output = nn.Linear(hidden_dim//2, 1)
        
    def forward(self, x):
        # Project to hidden dimension
        x = self.input_proj(x)  # [batch_size, hidden_dim]
        
        # Calculate attention weights
        attention_weights = torch.softmax(self.attention(x), dim=0)  # [batch_size, 1]
        
        # Apply attention
        attended_x = x * attention_weights  # [batch_size, hidden_dim]
        
        # Feedforward network
        ff_output = self.feed_forward(attended_x)
        
        # Output prediction
        return self.output(ff_output).squeeze(-1)

class AutoencoderClassifier(nn.Module):
    def __init__(self, input_dim, latent_dim=16, dropout_rate=0.3):
        super().__init__()
        # Encoder - smaller hidden layers
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 40),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(40, 28),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(28, latent_dim)
        )
        
        # Decoder - correspondingly reduced
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 28),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(28, 40),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(40, input_dim)
        )
        
        # Classifier - simplified
        self.classifier = nn.Sequential(
            nn.Linear(latent_dim, 8),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(8, 1)
        )
    
    def forward(self, x, return_reconstruction=False):
        encoded = self.encoder(x)
        
        if return_reconstruction:
            decoded = self.decoder(encoded)
            classification = self.classifier(encoded).squeeze(-1)
            return classification, decoded
        else:
            return self.classifier(encoded).squeeze(-1)

# Function to count model parameters
def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def print_model_sizes(input_dim=50):
    """Print parameter count comparison of three models"""
    print(f"\nModel Parameter Count Comparison (input_dim={input_dim}):")
    print("="*50)
    
    models = {
        'SimpleDNN': SimpleDNN(input_dim),
        'AttentionMLP': AttentionMLP(input_dim), 
        'AutoencoderClassifier': AutoencoderClassifier(input_dim)
    }
    
    for name, model in models.items():
        param_count = count_parameters(model)
        print(f"{name:20}: {param_count:,} parameters")
    
    print("="*50 + "\n")

# 4. Training and evaluation functions
def train_model(model, loader, criterion, optimizer, device, is_autoencoder=False):
    model.train()
    for x, y in loader:
        x, y = x.to(device), y.to(device)
        optimizer.zero_grad()
        
        if is_autoencoder and hasattr(model, 'decoder'):
            # AutoencoderClassifier training
            classification, reconstruction = model(x, return_reconstruction=True)
            class_loss = criterion(classification, y)
            recon_loss = nn.MSELoss()(reconstruction, x)
            loss = class_loss + 0.1 * recon_loss  # Balance classification and reconstruction loss
        else:
            # Standard model training
            out = model(x)
            loss = criterion(out, y)
        
        loss.backward()
        optimizer.step()

def eval_model(model, loader, device):
    model.eval()
    y_true, y_pred, y_prob = [], [], []
    with torch.no_grad():
        for x, y in loader:
            x = x.to(device)
            out = model(x)
            prob = torch.sigmoid(out).cpu().numpy()
            pred = (prob > 0.5).astype(int)
            y_true.extend(y.cpu().numpy())
            y_pred.extend(pred)
            y_prob.extend(prob)
    return np.array(y_true), np.array(y_pred), np.array(y_prob)

# 5. Metrics calculation
def calc_metrics(y_true, y_pred, y_prob):
    auc = roc_auc_score(y_true, y_prob) if len(np.unique(y_true)) > 1 else np.nan
    acc = accuracy_score(y_true, y_pred)
    sen = recall_score(y_true, y_pred, pos_label=1)
    spc = recall_score(y_true, y_pred, pos_label=0)
    return {'auc': auc, 'acc': acc, 'sen': sen, 'spc': spc}

# 6. Single experiment function
def run_single_experiment(model_class, model_name, X_source, y_source, X_target, y_target, 
                         strategy, fold, train_idx, test_idx, device, epochs=50, batch_size=32):
    """Run single experiment"""
    is_autoencoder = model_name == 'AutoencoderClassifier'
    
    if strategy == 'Regional':
        # Train only on source domain, test on target domain
        model = model_class(X_source.shape[1]).to(device)
        optimizer = optim.Adam(model.parameters(), lr=1e-3)
        criterion = nn.BCEWithLogitsLoss()
        ds = TensorDataset(torch.tensor(X_source, dtype=torch.float32), torch.tensor(y_source, dtype=torch.float32))
        loader = DataLoader(ds, batch_size=batch_size, shuffle=True)
        for ep in range(epochs):
            train_model(model, loader, criterion, optimizer, device, is_autoencoder)
        # Testing
        test_ds = TensorDataset(torch.tensor(X_target[test_idx], dtype=torch.float32), torch.tensor(y_target[test_idx], dtype=torch.float32))
        test_loader = DataLoader(test_ds, batch_size=batch_size)
        y_true, y_pred, y_prob = eval_model(model, test_loader, device)
        
    elif strategy == 'Transfer':
        # First train on source domain, then fine-tune on target domain training set
        model = model_class(X_source.shape[1]).to(device)
        optimizer = optim.Adam(model.parameters(), lr=1e-3)
        criterion = nn.BCEWithLogitsLoss()
        ds = TensorDataset(torch.tensor(X_source, dtype=torch.float32), torch.tensor(y_source, dtype=torch.float32))
        loader = DataLoader(ds, batch_size=batch_size, shuffle=True)
        for ep in range(epochs):
            train_model(model, loader, criterion, optimizer, device, is_autoencoder)
        # Fine-tuning
        finetune_ds = TensorDataset(torch.tensor(X_target[train_idx], dtype=torch.float32), torch.tensor(y_target[train_idx], dtype=torch.float32))
        finetune_loader = DataLoader(finetune_ds, batch_size=batch_size, shuffle=True)
        for ep in range(epochs):
            train_model(model, finetune_loader, criterion, optimizer, device, is_autoencoder)
        # Testing
        test_ds = TensorDataset(torch.tensor(X_target[test_idx], dtype=torch.float32), torch.tensor(y_target[test_idx], dtype=torch.float32))
        test_loader = DataLoader(test_ds, batch_size=batch_size)
        y_true, y_pred, y_prob = eval_model(model, test_loader, device)
        
    elif strategy == 'Independent':
        # Train only on target domain training set
        model = model_class(X_target.shape[1]).to(device)
        optimizer = optim.Adam(model.parameters(), lr=1e-3)
        criterion = nn.BCEWithLogitsLoss()
        train_ds = TensorDataset(torch.tensor(X_target[train_idx], dtype=torch.float32), torch.tensor(y_target[train_idx], dtype=torch.float32))
        train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
        for ep in range(epochs):
            train_model(model, train_loader, criterion, optimizer, device, is_autoencoder)
        # Testing
        test_ds = TensorDataset(torch.tensor(X_target[test_idx], dtype=torch.float32), torch.tensor(y_target[test_idx], dtype=torch.float32))
        test_loader = DataLoader(test_ds, batch_size=batch_size)
        y_true, y_pred, y_prob = eval_model(model, test_loader, device)
    
    metrics = calc_metrics(y_true, y_pred, y_prob)
    return metrics['auc']

# 7. Permutation Test function
def permutation_test(model_class, model_name, X_source, y_source, X_target, y_target,
                    strategy, device, n_permutations=50, n_repeats=5, n_folds=3, epochs=50):
    """
    Execute permutation test
    Returns real AUC and p-value
    """
    print(f"  Running permutation test for {model_name} - {strategy}...")
    
    # 1. Calculate real performance
    real_aucs = []
    for repeat in range(n_repeats):
        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=repeat)
        for fold, (train_idx, test_idx) in enumerate(skf.split(X_target, y_target)):
            auc = run_single_experiment(model_class, model_name, X_source, y_source, 
                                      X_target, y_target, strategy, fold, train_idx, test_idx, device, epochs)
            if not np.isnan(auc):
                real_aucs.append(auc)
    
    real_auc_mean = np.mean(real_aucs) if real_aucs else np.nan
    
    # 2. Permutation testing
    permuted_aucs = []
    for perm in tqdm(range(n_permutations), desc=f"    Permutation", leave=False):
        # Shuffle labels
        y_target_permuted = np.random.permutation(y_target)
        
        perm_aucs = []
        for repeat in range(n_repeats):
            skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=repeat)
            for fold, (train_idx, test_idx) in enumerate(skf.split(X_target, y_target_permuted)):
                try:
                    auc = run_single_experiment(model_class, model_name, X_source, y_source,
                                              X_target, y_target_permuted, strategy, fold, 
                                              train_idx, test_idx, device, epochs)
                    if not np.isnan(auc):
                        perm_aucs.append(auc)
                except:
                    continue
        
        if perm_aucs:
            permuted_aucs.append(np.mean(perm_aucs))
    
    # 3. Calculate p-value
    if permuted_aucs and not np.isnan(real_auc_mean):
        p_value = np.mean([perm_auc >= real_auc_mean for perm_auc in permuted_aucs])
    else:
        p_value = np.nan
    
    return real_auc_mean, p_value, len(permuted_aucs)

# 8. Main loop
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# Model selection configuration
MODEL_TYPES = {
    'SimpleDNN': SimpleDNN,
    'AttentionMLP': AttentionMLP,
    'AutoencoderClassifier': AutoencoderClassifier
}

# Select models to test
models_to_test = ['SimpleDNN', 'AttentionMLP', 'AutoencoderClassifier']

# Experiment parameters
n_permutations = 50
n_repeats = 5  # Reduce repetitions to speed up
n_folds = 3
epochs = 30  # Reduce training epochs to speed up

# Print model parameter count comparison
print_model_sizes(input_dim=50)

# Store permutation test results
permutation_results = []

print("Starting Permutation Tests...")
print("="*80)

for source_host in tqdm(host_map, desc="Source host"):
    source_keys = host_map[source_host]
    X_source = np.concatenate([profiles[k].values for k in source_keys])
    y_source = np.concatenate([labels[k] for k in source_keys])
    y_source = (y_source == 'sick').astype(int)
    target_hosts = [h for h in host_map if h != source_host]
    
    for target_host in tqdm(target_hosts, desc=f"Target host ({source_host})", leave=False):
        for target_key in tqdm(host_map[target_host], desc=f"Target dataset ({target_host})", leave=False):
            X_target = profiles[target_key].values
            y_target = (labels[target_key] == 'sick').astype(int)
            
            # Check data quality
            if len(np.unique(y_target)) < 2:
                print(f"Skipping {target_key}: insufficient class diversity")
                continue
            
            for model_name in tqdm(models_to_test, desc=f"Model ({target_key})", leave=False):
                model_class = MODEL_TYPES[model_name]
                
                for strategy in ['Regional', 'Transfer', 'Independent']:
                    # Execute permutation test
                    real_auc, p_value, n_valid_perms = permutation_test(
                        model_class, model_name, X_source, y_source, X_target, y_target,
                        strategy, device, n_permutations, n_repeats, n_folds, epochs
                    )
                    
                    # Record results
                    result = {
                        'source_host': source_host,
                        'target_host': target_host,
                        'target_dataset': target_key,
                        'model': model_name,
                        'strategy': strategy,
                        'real_auc': real_auc,
                        'p_value': p_value,
                        'n_permutations': n_permutations,
                        'n_valid_permutations': n_valid_perms,
                        'significant': p_value < 0.05 if not np.isnan(p_value) else False
                    }
                    permutation_results.append(result)
                    
                    # Real-time output results
                    sig_marker = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
                    print(f"    {model_name} - {strategy}: AUC={real_auc:.4f}, p={p_value:.4f} {sig_marker}")

# 9. Save results
results_df = pd.DataFrame(permutation_results)
results_df.to_csv('permutation.csv', index=False)

# 10. Generate summary report
print("\n" + "="*80)
print("PERMUTATION TEST SUMMARY")
print("="*80)

# Significance results statistics
significant_results = results_df[results_df['significant'] == True]
print(f"\nSignificant results (p < 0.05): {len(significant_results)} / {len(results_df)}")

# Significance statistics by model
model_significance = results_df.groupby('model').agg({
    'significant': ['sum', 'count'],
    'p_value': 'mean',
    'real_auc': 'mean'
}).round(4)

print("\nModel Performance Summary:")
print(model_significance)

# Significance statistics by strategy
strategy_significance = results_df.groupby('strategy').agg({
    'significant': ['sum', 'count'],
    'p_value': 'mean',
    'real_auc': 'mean'
}).round(4)

print("\nStrategy Performance Summary:")
print(strategy_significance)

# Best results
best_results = results_df.nlargest(10, 'real_auc')[['model', 'strategy', 'source_host', 'target_host', 'real_auc', 'p_value', 'significant']]
print("\nTop 10 Results by AUC:")
print(best_results.round(4))

print(f"\nDetailed results saved to: permutation.csv")
print("="*80)